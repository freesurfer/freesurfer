/*============================================================================
  Copyright (c) 1996 Martin Sereno and Anders Dale
  ===========================================================================*/

// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2003/05/07 17:35:55 $
// Revision       : $Revision: 1.145 $
char *VERSION = "$Revision: 1.145 $";

#define TCL
#define TKMEDIT 
/*#if defined(Linux) || defined(sun) || defined(SunOS) */
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
/*#endif */
#include <fcntl.h>
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "const.h"
#include "transform.h"
#include "version.h"

#include "tkmedit.h"


#define SET_TCL_ENV_VAR 1
#include <tcl.h>
// #include <tclDecls.h>
#include <tk.h>
#include <tix.h>
#include <blt.h>

/* Blt_Init is not declared in the blt.h for some reason, so this is just
   to keep the compiler for bitching. */
#ifndef Blt_Init
int Blt_Init ( Tcl_Interp* interp );
#endif 
#ifndef Blt_SafeInit
int Blt_SafeInit ( Tcl_Interp* interp );
#endif 
#ifndef Tix_Init
int Tix_Init ( Tcl_Interp* interp );
#endif 
#ifndef Tix_SafeInit
int Tix_SafeInit ( Tcl_Interp* interp );
#endif 

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "MRIio_old.h"
#include "volume_io.h"
#include "rgb_image.h"
#include "fio.h"
#include "mri_conform.h"

#ifndef OPENGL
#define OPENGL
#endif

#include "proto.h"
#include "macros.h"
#include "matrix.h"

MATRIX *gm_screen2ras = NULL ;
MATRIX *gm_ras2screen = NULL ;

char *tkm_ksaErrorStrings [tkm_knNumErrorCodes] = {
  "No error.",
  "A parameter to this function was invalid.",
  "The environment variable MRI_DIR is not defined.",
  "The environment variable SUBJECTS_DIR is not defined.",
  "Couldn't get the current working directory.",
  "Couldn't initialize Tcl.",
  "Couldn't initialize Tk.",
  "Couldn't initialize Tix.",
  "Couldn't initialize BLT.",
  "Couldn't read the anatomical volume.",
  "Couldn't load the surface.",
  "Couldn't load the label.",
  "Couldn't load the surface vertex set.",
  "Couldn't load the color table.",
  "Couldn't load the head points list.",
  "Couldn't import the segmentation volume.",
  "Couldn't load the functional overlay.",
  "Couldn't load the functional time course.",
  "Couldn't load the transform.",
  "Couldn't load the GCA volume.",
  "Couldn't load the VLI volume.",
  "Couldn't load the DTI volume.",
  "Couldn't load the vector field.",
  "Error accessing a file.",
  "Error accessing the anatomical volume.",
  "Error accessing the segmentation.",
  "Error accessing the functional volume.",
  "Error accessing the list.",
  "Error accessing surface.",
  "Couldnt write a file.",
  "A memory allocation failed.",
  "Tried to call a function on an unloaded surface.",
  "Functional overlay is not loaded.",
  "GCA volume is not loaded.",
  "Segmentation data is not loaded.",
  "The FA and EV volumes are different sizes.",
  "Couldn't cache the script name.",
  "gettimeofday failed.",
  "Unrecoverable error from an inner function."
};



#define NUM_UNDOS   257

#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))


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


static int zf, ozf;
static float fsf;

static int editflag = TRUE;
static int surflinewidth = 1;
// int editedimage = FALSE;


#include "xDebug.h"
#include "xTypes.h"
#include "xUtilities.h"
#include "xVoxel.h"
#include "xList.h"
#include "x3DList.h"
#include "tkmMeditWindow.h"


// ============================================================= FILENAME MGMT

static char *gsCommandLineSubjectsDir = NULL ;
/* is $SUBJECTS_DIR/subject_name if they used the "subject image" argument, 
   or the -f arg if they used that, or the cwd */
static char gsSubjectHomeDir[tkm_knPathLen] = ""; 
static char gsUserHomeDir[tkm_knPathLen]    = ""; /* cwd from shell */

static char gsTkTimerFileName[tkm_knPathLen] = "tktimer.data";

tkm_tErr SetSubjectHomeDirFromEnv ( char* isSubject );
tkm_tErr SetSubjectHomeDir    ( char* isHomeDir );
tkm_tErr FindUserHomeDir    ();

typedef enum {
  tkm_tFileName_Functional = 0,
  tkm_tFileName_Segmentation,
  tkm_tFileName_HeadPoints,
  tkm_tFileName_Surface,
  tkm_tFileName_Volume,
  tkm_tFileName_VolumeTransform,
  tkm_tFileName_Label,
  tkm_tFileName_GCA,
  tkm_tFileName_VLI,
  tkm_tFileName_RGB,
  tkm_tFileName_ControlPoints,
  tkm_tFileName_Edit,
  tkm_tFileName_TclScript,
  tkm_knNumFileNameTypes
} tkm_tFileName;

/* subdirectories local to subject's home dir */
char *ksaFileNameSubDirs[tkm_knNumFileNameTypes] = {
  "fmri", "", "bem", "surf", "mri", "mri/transforms", "label", "mri", "", "image/rgb", "tmp", "tmp", "lib/tcl"
};

/* input starts with ., gsUserHomeDir will be prepended. if ~ or nothing,
   gsSubjectHomeDir will be prepended. anything else, first use the subject
   home dir and then the subdir of the proper type, then the proper
   remainder of the input file name. if it starts with a /, nothing will be
   added. */
void MakeFileName ( char*     isInput,
		    tkm_tFileName  iType, 
		    char*         osCompleteFileName,
		    int         inDestSize );

/* attempts to extract a name from the data */
void ExtractSubjectName ( char* isDataSource,
			  char* osSubjectName );
void ExtractVolumeName  ( char* isDataSource,
			  char* osVolumeName );

/* only tries to prepend subdirectories if this flag is true. it is set to
   false if -f is used on the command line. */
tBoolean gEnableFileNameGuessing = TRUE;


// ==========================================================================

// ================================================== SELECTING CONTROL POINTS

/* returns distance to nearest control point on the same plane. returns
   0 if there isn't one. */
float FindNearestControlPoint ( xVoxelRef   iAnaIdx, 
				mri_tOrientation iPlane,
				xVoxelRef*   opCtrlPt );

/* makes a copy of the ctrl pt and puts it in the ctrl pt list */
void NewControlPoint         ( xVoxelRef iCtrlPt,
			       tBoolean  ibWriteToFile );

/* removes ctrl pt from the list and deletes it */
void DeleteControlPoint         ( xVoxelRef ipCtrlPt );

/* reas the control.dat file, transforms all pts from RAS space 
   to voxel space, and adds them as control pts */
void ProcessControlPointFile ();
tBoolean gbParsedControlPointFile = FALSE;

/* writes all control points to the control.dat file in RAS space */
void WriteControlPointFile   ();

x3DListRef gControlPointList = NULL;

/* function for comparing voxels */
xList_tCompare CompareVoxels ( void* inVoxelA, void* inVoxelB );

// ===========================================================================

// ================================================================== SURFACES

#include "mriSurface.h"

mriSurfaceRef gSurface[tkm_knNumSurfaceTypes];

tkm_tErr LoadSurface          ( tkm_tSurfaceType iType, 
				char*     isName );
tkm_tErr LoadSurfaceVertexSet ( tkm_tSurfaceType iType, 
				Surf_tVertexSet   iSet,
				char*     fname );
void   UnloadSurface          ( tkm_tSurfaceType iType );

void   WriteSurfaceValues     ( tkm_tSurfaceType iType, 
				char*     isFileName );

void GotoSurfaceVertex                    ( Surf_tVertexSet iSurface, 
					    int inVertex );
void FindNearestSurfaceVertex             ( Surf_tVertexSet iSet );
void FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet iSet );
void AverageSurfaceVertexPositions        ( int inNumAverages );

// ===========================================================================

// ========================================================= SELECTING REGIONS

/* selecting regions works much like editing. the user chooses a brush
   shape and size and paints in a region. the tool can be toggled
   between selecting and unselecting. the selected pixels are kept in
   a voxel space for optimized retreival in the draw loop. there are
   the usual functions for adding and removing voxels as well as
   saving them out to a file. */

x3DListRef gSelectedVoxels;

tkm_tErr InitSelectionModule ();
void   DeleteSelectionModule ();

/* grabs the list of voxels selected and draws them into the buffer. */
void DrawSelectedVoxels ( char * inBuffer, int inPlane, int inPlaneNum );

/* adds or removes voxels to selections. if a voxel that isn't in the 
   selection is told to be removed, no errors occur. this is called from the 
   brush function. The iaAnaIdx parameter can be an array. */
void AddVoxelsToSelection      ( xVoxelRef  iaAnaIdx, int inCount );
void RemoveVoxelsFromSelection ( xVoxelRef  iaAnaIdx, int inCount );

/* clears the current selection */
void ClearSelection ();

/* write to and read from label files */
void SaveSelectionToLabelFile ( char * inFileName );
tkm_tErr LoadSelectionFromLabelFile ( char * inFileName );

/* send the selected voxels to the functional display module */
void GraphSelectedRegion ();

/* tell the functional volume to select all voxels contiguous to the cursor
   and having values above the func value at the cursor */
void SelectVoxelsByFuncValue ( FunV_tFindStatsComp iCompare );

/* values used in calculating selection color */
#define kSel_IntensifyValue 100
#define kSel_TooClose      10

// ===========================================================================

// ============================================================= VOLUME ACCESS

#include "mriVolume.h"

/* we declare tkm_knNumVolumeTypes, but we only use the main and aux */
static mriVolumeRef     gAnatomicalVolume[tkm_knNumVolumeTypes];
static int              gnAnatomicalDimensionX = 0;
static int              gnAnatomicalDimensionY = 0;
static int              gnAnatomicalDimensionZ = 0;
static mriTransformRef  gIdxToRASTransform    = NULL;
static tBoolean         gbAnatomicalVolumeDirty[tkm_knNumVolumeTypes];

tkm_tErr LoadVolume    ( tkm_tVolumeType iType,
			 char*     isFileName );

tkm_tErr UnloadVolume  ( tkm_tVolumeType iType );

/* saves the anatomical volume. if isPath is null, saves over the original. */
tkm_tErr SaveVolume    ( tkm_tVolumeType iType,
			 char*     isPath );


/* tells volume to load or unload this display transform */
tkm_tErr LoadDisplayTransform  ( tkm_tVolumeType iVolume,
				 char*      isFileName );
tkm_tErr UnloadDisplayTransform ( tkm_tVolumeType iVolume );

/* operations on the main volume */
void SnapshotVolume ();
void RestoreVolumeFromSnapshot ();

tkm_tErr SetVolumeDirty ( tkm_tVolumeType iVolume, tBoolean ibDirty );
tkm_tErr IsVolumeDirty ( tkm_tVolumeType iVolume, tBoolean* obDirty );

void SetVolumeBrightnessAndContrast ( tkm_tVolumeType iVolume,
				      float        ifBrightness,
				      float        ifContrast );

void ThresholdVolume ( int    inLevel,
		       tBoolean      ibAbove,
		       int        inNewLevel );
void FlipVolume       ( mri_tOrientation iAxis );
void RotateVolume    ( mri_tOrientation iAxis,
		       float      ifDegrees );

void EditAnatomicalVolumeInRangeArray ( tkm_tVolumeType iVolume,
					xVoxelRef       ipaVoxel, 
					int             inCount,
					Volm_tValue     inLow,
					Volm_tValue     inHigh, 
					Volm_tValue     inNewValue );

void SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
				 int             iAnaX0,
				 int             iAnaX1,
				 int             iAnaY0,
				 int             iAnaY1,
				 int             iAnaZ0,
				 int             iAnaZ1,
				 float           iNewValue );

int EditAnatomicalVolume ( xVoxelRef iAnaIdx, int inValue );
int EditAuxAnatomicalVolume ( xVoxelRef iAnaIdx, int inValue );

void ConvertRASToAnaIdx ( xVoxelRef iRAS,
			  xVoxelRef oAnaIdx );

void ConvertAnaIdxToRAS ( xVoxelRef iAnaIdx,
			  xVoxelRef oRAS );

float gfaBrightness[tkm_knNumVolumeTypes]; 
float gfaContrast[tkm_knNumVolumeTypes]; 

void SetVolumeBrightnessAndContrast  ( tkm_tVolumeType iVolume,
				       float        ifBrightness,
				       float        ifContrast );
void SendBrightnessAndContrastUpdate ( tkm_tVolumeType iVolume );

// ========================================================= FUNCTIONAL VOLUME

#include "tkmFunctionalVolume.h"

tkmFunctionalVolumeRef gFunctionalVolume = NULL;

tkm_tErr LoadFunctionalOverlay    ( char* isPathAndStem, 
            char* isOffsetPath,
            char* isRegistration );
tkm_tErr LoadFunctionalTimeCourse ( char* isPathAndStem, 
            char* isOffsetPath,
            char* isRegistration );

tkm_tErr SmoothOverlayData ( float ifSigma );

void TranslateOverlayRegisgtration ( float ifDistance, char isDirection );
void RotateOverlayRegistration     ( float ifDegrees, char isDirection );
void ScaleOverlayRegisgtration     ( float ifDistance, char isDirection );

// ===========================================================================

// ============================================================== SEGMENTATION

#include "mriColorLookupTable.h"

#define kfDefaultSegmentationAlpha 0.3

/* gSegmentationVolume is the normal segmentation volume,
   gPreviousSegmentationVolume is a backup that is created before it
   is recomputed in RecomputeSegmentation, and
   gSegmentationChangedVolume contains flags at the values that were
   changed in this session. */
static mriVolumeRef  gSegmentationVolume[tkm_knNumSegTypes];
static mriVolumeRef  gPreviousSegmentationVolume[tkm_knNumSegTypes];
static mriVolumeRef  gSegmentationChangedVolume[tkm_knNumSegTypes];

static mriColorLookupTableRef gColorTable[tkm_knNumSegTypes];
static float         gfSegmentationAlpha   = kfDefaultSegmentationAlpha;

/* Creates a new segementation volume from the settings from an
   existing anatomical volume, basically just copying it and erasing
   the value contents. */
tkm_tErr NewSegmentationVolume  ( tkm_tSegType    iVolume,
				  tkm_tVolumeType iFromVolume,
				  char*           inColorFileName );

/* Loads a segmentation volume, can be a different dimension than the
   anatomical (i think) */
tkm_tErr LoadSegmentationVolume ( tkm_tSegType    iVolume,
				  char*           inVolumeDir,
				  char*           inColorFileName );

/* Saves a segmentation volume to COR format. */
void SaveSegmentationVolume     ( tkm_tSegType    iVolume,
				  char*           inVolumeDir );

/* Uses the ChangedVolume to save a volume with only those label
   values that have been editied in this session. Basically an AND of
   the SegmentationVolume and ChangedVolume. */
typedef struct {
  mriVolumeRef mSrcVolume;
  mriVolumeRef mDestVolume;
} tkm_tExportSegmentationParams, *tkm_tExportSegmentationParamsRef;
tkm_tErr ExportChangedSegmentationVolume ( tkm_tSegType    iVolume,
					   char*           inVolumeDir );

/* Creates an segmentation from a surface annotation file. Only the
   voxel values that lie on the surface will be filled. */
tkm_tErr ImportSurfaceAnnotationToSegmentation ( tkm_tSegType iVolume,
						 char*    inAnnotationFileName,
						 char*    inColorFileName );

tkm_tErr LoadSegmentationColorTable ( tkm_tSegType iVolume,
				      char*        inColorFileName );


/* Sets the display alpha for determining the opacity of the
   segmentation overlay. */
void SetSegmentationAlpha ( float ifAlpha );

/* Gets a color at an index in a seg volume. Given a base color, this
   function will blend the label color according to the opacity alpha
   and return the new color. */
void GetSegmentationColorAtVoxel ( tkm_tSegType iVolume,
				   xVoxelRef   inVoxel,
				   xColor3fRef iBaseColor,
				   xColor3fRef oColor );

/* Returns the name of the label and/or label index for a location. */
void GetSegLabel        ( tkm_tSegType iVolume,
			  xVoxelRef   ipVoxel, 
			  int*        onIndex,
			  char*       osLabel );

/* Adds all voxels in a label to the selection. */
tkm_tErr SelectSegLabel ( tkm_tSegType iVolume,
			  int          inIndex );

/* Graphs the average of all the voxels in a label. */
tkm_tErr GraphSegLabel  ( tkm_tSegType iVolume,
		     int          inIndex );

/* Calls GCAreclassifyUsingGibbsPriors to recalc the segmentation
   using changed values, finally updating the segemention values with
   the changes. If gDisplayIntermediateResults is on, will draw
   intermediate volumes as well. The gRecaluclatingSegemention is used
   during this function by the callback to update the proper
   segmentation. */
tkm_tSegType gRecaluclatingSegemention = tkm_tSegType_Main;
void RecomputeSegmentation ( tkm_tSegType iVolume );
static void RecomputeUpdateCallback( MRI* iMRIValues );

/* A call back for a visit function. If the two values are equal
   (there's a float in ipnTarget) then this voxel will be added to the
   selection. */
Volm_tVisitCommand AddSimilarVoxelToSelection ( xVoxelRef iAnaIdx,
						float     iValue,
						void*     ipnTarget );
/* A call back for a visit function. If the two values are equal
   (there's a float in ipnTarget) then this voxel will be added to the
   functional selection so it can later be graphed. */
Volm_tVisitCommand AddSimilarVoxelToGraohAvg  ( xVoxelRef iAnaIdx,
						float     iValue,
						void*     ipnTarget );

/* A call back for a visit function. When run on the changed volume by
   ExportChangedSegmentationVolume, for each value that is 1, sets the
   value in the dest volume to the value in the src volume. Uses
   tkm_tExportSegmentationParamsRef. */
Volm_tVisitCommand SetChangedSegmentationValue ( xVoxelRef iAnaIdx,
						 float     iValue,
						 void*     iData );

/* A callback to an undo entry. Sets the segmentation volume value at
   this index, as well as the changed volume index.  */
int  EditSegmentation ( tkm_tSegType iVolume,
			xVoxelRef    iAnaIdx,
			int          inIndex );

/* Calculates the volume of a contiguous label. */
void CalcSegLabelVolume        ( tkm_tSegType iVolume,
			    xVoxelRef    ipVoxel,
			    int*         onVolume );
void IterateCalcSegLabelVolume ( tkm_tSegType iVolume,
			    xVoxelRef    iAnaIdx, 
			    int          inIndex, 
			    tBoolean*    iVisited, 
			    int*         onVolume );

/* Use to set the segmentation value. Also adds to the undo
   list. TODO: Add it to the undo list, it doesn't currently. */
void SetSegmentationValue    ( tkm_tSegType iVolume,
			       xVoxelRef    iAnaIdx,
			       int          inIndex );
void SetSegmentationValues   ( tkm_tSegType iVolume,
			       xVoxelRef    iaAnaIdx,
			       int          inCount,
			       int          inIndex );

/* Flood fills a segmentation volume with various parameters. */
typedef struct {
  int      mnTargetValue;      /* the value to look for */
  int      mnFuzziness;        /* the range around the target value */
  int      mnDistance;        /* maximum distance to affect */
  xVoxel   mAnchor;        /* the anchor of the fill */
  tBoolean mb3D;          /* whether to go in all directions */
  mri_tOrientation mOrientation;      /* plane to go in if 2D */
  tkm_tVolumeTarget mSource;        /* what volume to get the target val */
  int      mnNewIndex;        /* new value to set in the parc */
  int      mnIterationDepth;   /* how many iteration levels currently */
  int      mnIterationLimit;   /* how deep to go */
  tBoolean mbIterationLimitReached; /* whether we were stopped by the limit */
} tkm_tParcFloodFillParams, *tkm_tParcFloodFillParamsRef;

void FloodFillSegmentation        ( tkm_tSegType                iVolume,
				    xVoxelRef                   iAnaIdx,
				    tkm_tParcFloodFillParamsRef iParams );
void IterateFloodFillSegmentation ( tkm_tSegType                iVolume, 
				    xVoxelRef                   iAnaIdx,
				    tkm_tParcFloodFillParamsRef iParams,
				    tBoolean*                   iVisited );

#define knMaxFloodFillIteration 10000

// ===========================================================================

// =============================================================== DTI VOLUMES

mriVolumeRef gDTIVolume;
tAxis gaDTIAxisForComponent[xColr_knNumComponents];
float gfDTIAlpha = 1.0;

tkm_tErr LoadDTIVolume ( char*              isNameEV,
                         char*              isNameFA,
			 tAxis              iRedAxis,
			 tAxis              iGreenAxis,
			 tAxis              iBlueAxis );

void GetDTIColorAtVoxel ( xVoxelRef        iAnaIdx,
			  mri_tOrientation iPlane,
			  xColor3fRef      iBaseColor,
			  xColor3fRef      oColor );

void SetDTIAlpha ( float ifAlpha );

// ===========================================================================

// ============================================================ VECTOR VOLUMES

mriVolumeRef gVectorField;

tkm_tErr LoadVectorField ( char* isName );

void GetVectorAtVoxel ( xVoxelRef iAnaIdx,
      xVoxelRef oDirection );

// ===========================================================================

// ================================================================= UNDO LIST

#include "xUndoList.h"

/* this is a pretty simple implementation
   of undo that only supports pixel 
   editing. when editing, the pixels that
   were changed are saved along with their
   previous values in a list, one list per
   editing click. */
xUndoListRef gUndoList = NULL;
tkm_tErr InitUndoList ();
void   DeleteUndoList ();
void   ClearUndoList ();

/* protoype for the kind of function that
   gets called for each undo entry. should
   return the old value. */
typedef int(*tUndoActionFunction) ( xVoxelRef, int );

/* when pixels are editied, they are added
   to the list. if the user hits undo, the
   entire list is drawn to the screen, 
   using the SetVoxelValue() function. at
   the same time, a new list is made
   to save the positions of all the restored
   voxel values. that list becomes the
   new undo list, so you can effectivly
   undo an undo. */
void AddVoxelAndValueToUndoList ( tUndoActionFunction iFunction,
          xVoxelRef        iVoxel, 
          int          inValue );
void RestoreUndoList ();

/* we need a struct for the undo list. this
   is what we add to it and what we get
   back when the list is restored. */
typedef struct {
  xVoxelRef        mVoxel;   /* which voxel to set */
  int          mnValue;   /* the value to set it to */
  tUndoActionFunction mFunction; /* the function to call to do it */
} UndoEntry, *UndoEntryRef;

void NewUndoEntry      ( UndoEntryRef*    opEntry, 
            tUndoActionFunction iFunction,
            xVoxelRef    iVoxel, 
            int      inValue );
void DeleteUndoEntry      ( UndoEntryRef*    ioEntry );
void PrintUndoEntry      ( UndoEntryRef    iEntry );

/* these are our callback functions for the
   undo list. the first deletes an entry
   and the second actually performs the
   undo action and hands back the undone
   voxel. */
void DeleteUndoEntryWrapper ( xUndL_tEntryPtr* ipEntryToDelete );
void UndoActionWrapper      ( xUndL_tEntryPtr  iUndoneEntry, 
            xUndL_tEntryPtr* opNewEntry );
void PrintEntryWrapper      ( xUndL_tEntryPtr  iEntry );


// ==========================================================================

/* =========================================================== UNDO VOLUME */

#include "xSparseVolume.h"

/* a volume for saving multiple edited voxel values so that they can be undone
   in a flood-fill style */

xSparseVolumeRef gUndoVolume = NULL;

typedef struct {
  int mRestoreValue;
} UndoVolumeEntry, *UndoVolumeEntryRef;

UndoVolumeEntryRef NewUndoVolumeEntry    ( int         iValue );
void       DeleteUndoVolumeEntry  ( UndoVolumeEntryRef*iEntry );
void       DeleteUndoVolumeEntryWrapper ( void*         iEntry );

/* initializes and deletes the volume */
tkm_tErr InitUndoVolume    ();
void   DeleteUndoVolume ();

/* adds a value to the volume at an anatomical index */
void   AddAnaIdxAndValueToUndoVolume ( xVoxelRef iAnaIdx,
           int       iValue );

/* sees if there is a value for this ana idx, i.e. if it can be undone */
tBoolean IsAnaIdxInUndoVolume         ( xVoxelRef iAnaIdx );

/* resotres the values for all voxels touching this one */
void   RestoreUndoVolumeAroundAnaIdx ( xVoxelRef iAnaIdx );

/* clears the volume */
void   ClearUndoVolume ();

/* ======================================================================= */

/* =========================================================== HEAD POINTS */

#include "mriHeadPointList.h"

mriHeadPointListRef gHeadPoints = NULL;

tkm_tErr LoadHeadPts ( char* isHeadPtsFile, 
           char* isTransformFile );

void RestoreHeadPoints        ();
void WriteHeadPointsTransform ();
void WriteHeadPointsFile      ();
void RotateHeadPts        ( float ifDegrees, 
        char  isDirection );
void TranslateHeadPts        ( float ifDistance, 
        char  isDirection );

void SetSelectedHeadPointLabel ( char* isNewLabel );

void AlignSelectedHeadPointToAnaIdx ( xVoxelRef iAnaIdx );

/* ======================================================================= */

/* =================================================================== GCA */

#include "gca.h"

GCA* gGCAVolume     = NULL;
TRANSFORM* gGCATransform = NULL;
int gDisplayIntermediateResults = True ;

tkm_tErr SaveGCA   ( char* isGCAFileName);
tkm_tErr LoadGCA   ( char* isGCAFileName, char* isTransformFileName );
tkm_tErr LoadGCARenormalization   ( char* isRenormFname );
tkm_tErr DeleteGCA ();

/* ======================================================================= */

/* =================================================================== VLI */

#include "vlabels.h"

VLI* gVLI1 = NULL;
VLI* gVLI2 = NULL;

tkm_tErr LoadVLIs ( char* isFileName1, char* isFileName2 );
tkm_tErr DeleteVLIs ();

/* ======================================================================= */

/* ========================================================== MEDIT WINDOW */

tkmMeditWindowRef gMeditWindow = NULL;

/* ======================================================================= */

/* ================================================================= TIMER */

#define ksTkTimerDataFileEnvVar   "TKTIMER_FILENAME"
#define ksTkTimerAutoStartEnvVar "TKTIMER_AUTOSTART"

static long  gnStartTime = 0;
static tBoolean gbTimerOn   = FALSE;

tkm_tErr StartTimer ();
tkm_tErr StopTimer  ();

/* ======================================================================= */

// ====================================================================== MISC


/* goto this subject's saved cursor
   in edit.dat */
void GotoSavedCursor ();

/* reads a saved cursor from an edit.dat
   file and set the current cursor to
   it */
void ReadCursorFromEditFile ( char* isFileName );

/* segfault handler. saves data to /temp. */
void HandleSegfault ( int );

/* where to find the interface script */
char gInterfaceScriptName [tkm_knPathLen] = "";

/* set and get the tcl interp to send
   the msg to */
void SetTclInterp ( Tcl_Interp * inInterp );
Tcl_Interp * GetTclInterp ();

/* send a tcl command */
void SendTCLCommand ( char * inCommand );


#define tkm_knMaxNumCachedCommands 500
char gCachedTclCommands[tkm_knMaxNumCachedCommands][tkm_knPathLen];
int gNumCachedCommands = 0;
void SendCachedTclCommands ();
void PrintCachedTclErrorDlogsToShell ();
/* the tcl interpreter */
Tcl_Interp * gTclInterp = NULL;

/* don't start accepting tcl commands
   until setup is complete */
static tBoolean gbAcceptingTclCommands = FALSE;

static tBoolean gbUseCsurfInterface = FALSE;

tkm_tErr EnqueueScript ( char* isFileName );
tkm_tErr ExecuteQueuedScripts ();

#define tkm_knMaxNumCachedScripts 10
char gsCachedScriptNames[tkm_knMaxNumCachedScripts][tkm_knPathLen];
int gnNumCachedScripts = 0;

// ===========================================================================


#ifdef Linux
extern void scale2x(int, int, unsigned char *);
#endif


void ParseCmdLineArgs( int argc, char *argv[] );
void WriteVoxelToControlFile ( xVoxelRef inVolumeVox );
void WriteVoxelToEditFile    ( xVoxelRef inVolumeVox );

void rotate_brain(float a,char c) ;
void translate_brain(float a,char c) ;
void UpdateAndRedraw ();
void pix_to_rgb(char *fname) ; // another method of saving a screenshot
void scrsave_to_rgb(char *fname) ; // use scrsave to save a screen shot
void save_rgb(char *fname) ; // saves a screen shot
//void edit_pixel(int action);
void write_dipoles(char *fname) ;
void write_decimation(char *fname) ;
void read_hpts(char *fname) ;
void read_htrans(char *fname) ;
void write_htrans(char *fname) ;
void smooth_3d(int niter) ;
void wmfilter_corslice(int imc) ;
void alloc_second_im(void) ;
void smooth_surface(int niter) ;

char *Progname ;

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


void ParseCmdLineArgs ( int argc, char *argv[] ) {
  
  tkm_tErr      eResult          = tkm_tErr_NoErr;
  FunV_tErr      eFunctional          = FunV_tErr_NoError;
  int        nNumProcessedVersionArgs = 0;
  int        nCurrentArg        = 0;
  int        bFatalError        = FALSE;
  char        sArg[tkm_knPathLen]      = "";
  char        sError[tkm_knErrStringLen]    = "";
  tBoolean      bSubjectDeclared      = FALSE;
  tBoolean      bUsingMRIRead        = FALSE;
  char        sSubject[tkm_knPathLen]    = "";
  char        sImageDir[tkm_knPathLen]    = "";
  tBoolean      bSurfaceDeclared      = FALSE;
  char        sSurface[tkm_knPathLen]    = "";
  tBoolean      bLocalImageDir      = FALSE;
  tBoolean      bNoEdit        = FALSE;
  tBoolean      bLoadingAuxVolume      = FALSE;
  char        sAuxVolume[tkm_knPathLen]    = "";
  tBoolean      bLoadingMainTransform      = FALSE;
  char        sMainTransform[tkm_knPathLen]    = "";
  tBoolean      bLoadingAuxTransform      = FALSE;
  char        sAuxTransform[tkm_knPathLen]    = "";
  tBoolean      bLoadingOverlay      = FALSE;
  char        sOverlayPathAndStem[tkm_knPathLen]  = "";
  char        sOverlayOffsetPathAndStem[tkm_knPathLen] = "";
  char        sOverlayRegistration[tkm_knPathLen]  = "";
  char*        psOverlayOffsetPathAndStem    = NULL;
  char*        psOverlayRegistration      = NULL;
  tBoolean      bLoadingTimeCourse      = FALSE;
  char        sTimeCoursePathAndStem[tkm_knPathLen] = "";
  char        sTimeCourseOffsetPathAndStem[tkm_knPathLen]  = "";
  char        sTimeCourseRegistration[tkm_knPathLen]= "";
  char*        psTimeCourseOffsetPathAndStem    = NULL;
  char*        psTimeCourseRegistration    = NULL;
  tBoolean      bEnablingRegistration      = FALSE;
  tBoolean      bLoadingSegmentation      = FALSE;
  char        sSegmentationPath[tkm_knPathLen]  = "";
  char        sSegmentationColorFile[tkm_knPathLen] = "";
  tBoolean      bLoadingAuxSegmentation      = FALSE;
  char        sAuxSegmentationPath[tkm_knPathLen]  = "";
  char        sAuxSegmentationColorFile[tkm_knPathLen] = "";
  tBoolean      bLoadingVLI       = FALSE;
  char          sVLIFile1[tkm_knPathLen] = "";
  char          sVLIFile2[tkm_knPathLen] = "";
  tBoolean      bThresh        = FALSE;
  FunV_tFunctionalValue    min          = 0;
  tBoolean      bMid          = FALSE;
  FunV_tFunctionalValue    mid          = 0;
  tBoolean      bSlope        = FALSE;
  FunV_tFunctionalValue    slope          = 0;
  tBoolean      bSmooth        = FALSE;
  float        smoothSigma        = 0;
  tBoolean      bRevPhaseFlag        = FALSE;
  int        nRevPhaseFlag        = 0;
  tBoolean      bTruncPhaseFlag      = FALSE;
  int        nTruncPhaseFlag      = 0;
  tBoolean      bUseOverlayCacheFlag      = FALSE;
  int        nUseOverlayCacheFlag      = 0;
#if 0
  tBoolean      bSetConversionMethod      = FALSE;
#endif
#if 0
  FunD_tConversionMethod  convMethod     = FunD_tConversionMethod_FFF;
#endif
  tBoolean      bLoadingHeadPts      = FALSE;
  tBoolean      bHaveHeadPtsTransform      = FALSE;
  char        sHeadPts[tkm_knPathLen]    = "";
  char        sHeadPtsTransform[tkm_knPathLen]  = "";
  char        sLabel[tkm_knPathLen]      = "";
  tBoolean      bLoadingLabel        = FALSE;
  char        sSubjectTest[tkm_knPathLen]    = "";
  char        sScriptName[tkm_knPathLen]    = "";
  float        fSegmentationAlpha      = 0;
  tBoolean      bSegmentationAlpha      = FALSE;
  float         fBrightnessMain        = 0;
  float         fContrastMain          = 0;
  tBoolean      bBrightContrastMain    = FALSE;
  float         fBrightnessAux        = 0;
  float         fContrastAux          = 0;
  tBoolean      bBrightContrastAux    = FALSE;

  DebugEnterFunction( ("ParseCmdLineArgs( argc=%d, argv=%s )", 
           argc, argv[0]) );
  
  /* first get the functional threshold so we don't overwrite the defaults */
  DebugNote( ("Getting default functional threshold") );
  eFunctional = FunV_GetThreshold( gFunctionalVolume, &min, &min, &slope );
  DebugAssertThrow( (FunV_tErr_NoError == eFunctional ) );
  
  /* First look for the version option and handle that. If found,
     shorten our argc and argv count. If those are the only args we
     had, exit. */
  /* rkt: check for and handle version tag */
  nNumProcessedVersionArgs = handle_version_option (argc, argv, "$Id: tkmedit.c,v 1.145 2003/05/07 17:35:55 kteich Exp $");
  if (nNumProcessedVersionArgs && argc - nNumProcessedVersionArgs == 1)
    exit (0);
  argc -= nNumProcessedVersionArgs;

  if (argc<2) {
    
    printf("usage: tkmedit {[subject image_type]|[-f absolute_path]}\n");
    printf("         [surface_file] \n");
    printf("         [options ...]\n");
    printf("\n");
    printf("Anatomical Data\n");
    printf("\n");
    printf("subject image_type  : looks in $SUBJECTS_DIR/subject/mri for image_type\n");
    printf("-f absolute_path    : specify volume directory or file\n");
    printf("\n");
    printf("Surface\n");
    printf("\n");
    printf("surface_file   : surface file to load (relative to $SUBJECTS_DIR/surf \n");
    printf("               : or absolute)\n");
    printf("\n");
    printf("Options\n");
    printf("\n");
    printf("-aux <volume>  : load volume as auxilliary anatomical volume. relative to\n");
    printf("               : in $SUBJECTS_DIR/subject/mri or specify absolute path\n");
    printf("\n");
    printf("-bc-main <brightness> <contrast> : brightness and contrast for main volume\n");
    printf("-bc-aux <brightness> <contrast> : brightness and contrast for aux volume\n");
    printf("\n");
    printf("-overlay <path/stem>      : load functional overlay volume\n");
    printf("-overlay-reg <registration> : load registration file for overlay volume \n");
    printf("                            : (default is register.dat in same path as\n");
    printf("                            :  volume)\n");
    printf("\n");
    printf("-fthresh <value>      : specfify min, mid, and slope threshold\n");
    printf("-fmid <value>        : values for functional overlay display\n");
    printf("-fslope <value>        : (default is 0, 1.0, and 1.0)\n");
    printf("-fsmooth <sigma>      : smooth functional overlay after loading\n");
    printf("\n");
    printf("-revphaseflag <1|0>      : reverses phase display in overlay (default off)\n");
    printf("-truncphaseflag <1|0>      : truncates overlay values below 0 (default off)\n");
    printf("-overlaycache <1|0>      : uses overlay cache (default off)\n");
    printf("\n");
    printf("-sdir <subjects dir>       : (default is getenv(SUBJECTS_DIR)\n");
    printf("-timecourse <path/stem>    : load functional timecourse volume\n");
    printf("-timecourse-reg <registration>  : load registration file for timecourse   \n");
    printf("                                : volume (default is register.dat in\n");
    printf("                                : same path as volume)\n");
    printf("-timecourse-offset <path/stem>  : load timecourse offset volume\n");
    printf("\n");
    printf("-segmentation <volume> <colors>    : load segmentation volume and color file\n");
    printf("-aux-segmentation <volume> <colors>: load aux segmentation volume and color file\n");
    printf("-segmentation-opacity <opacity>    : opacity of the segmentation \n");
    printf("                                  : overlay (default is 0.3)\n");
    printf("\n");
    printf("-headpts <points> [<trans>]   : load head points file and optional\n");
    printf("                              : transformation\n");
    printf("\n");
    printf("-interface script    : scecify interface script (default is tkmedit.tcl)\n");
    printf("\n");
    printf("--version            : print version of tkmedit\n");
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
    DebugNote( ("Copying argv[%d]", nCurrentArg) );
    xUtil_strncpy( sArg, argv[nCurrentArg], sizeof(sArg) );
    
    /* check for a option */
    if( '-' == sArg[0] ) {
      
      if ( MATCH( sArg, "-tcl" ) ) {
  
	/* check for one more. */
	if( argc > nCurrentArg + 1 
	    && '-' != argv[nCurrentArg+1][0] ) {
    
	  /* copy the script name in */
	  DebugNote( ("Parsing -tcl option") );
	  xUtil_strncpy( sScriptName, argv[nCurrentArg+1], 
			 sizeof(sScriptName) );
    
	  /* queue the script */
	  EnqueueScript( sScriptName );
	  
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -tcl option",
			    "Expected an argument",
			    "This option needs an argument: the name of the "
			    "script to run." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-bc-main" ) ) {
  
	/* check for the 2 values following the switch */
	if( argc > nCurrentArg + 2
	    && '-' != argv[nCurrentArg+1][0]
	    && '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -bc-main option") );
	  fBrightnessMain = atof( argv[nCurrentArg+1] );
	  fContrastMain = atof( argv[nCurrentArg+2] );
	  bBrightContrastMain = TRUE;
	  nCurrentArg +=3 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -bc-main option",
			    "Expected two arguments",
			    "This option needs two arguments: the brightness"
			    "and the contrast for the volume." );
	  nCurrentArg += 1;
	}
	
      } else if( MATCH( sArg, "-bc-aux" ) ) {
  
	/* check for the 2 values following the switch */
	if( argc > nCurrentArg + 2
	    && '-' != argv[nCurrentArg+1][0]
	    && '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -bc-aux option") );
	  fBrightnessAux = atof( argv[nCurrentArg+1] );
	  fContrastAux = atof( argv[nCurrentArg+2] );
	  bBrightContrastAux = TRUE;
	  nCurrentArg +=3 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -bc-aux option",
			    "Expected two arguments",
			    "This option needs two arguments: the brightness"
			    "and the contrast for the volume." );
	  nCurrentArg += 1;
	}
	
      } else if( MATCH( sArg, "-o" ) ) {
	
	/* make sure there are enough args */
	if( argc > nCurrentArg + 2 &&
	    '-' != argv[nCurrentArg+1][0] &&
	    '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* read the overlay path and stem */
	  DebugNote( ("Parsing overlay in -o option") );
	  xUtil_snprintf( sOverlayPathAndStem, sizeof(sOverlayPathAndStem),
			  "%s/%s",argv[nCurrentArg+1], argv[nCurrentArg+2] );
	  bLoadingOverlay = TRUE;
	  nCurrentArg += 3;
	  
	} else {
	  
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -o option",
			    "Expected an argument",
			    "This option needs two arguments: the path and "
			    "the stem of the data." );
	  nCurrentArg ++;
	}
  
	/* check for two more. */
	if( argc > nCurrentArg + 1 
	    && '-' != argv[nCurrentArg][0] ) {
	  
	  /* read in time course path and stem. */
	  DebugNote( ("Parsing time course in -o option") );
	  xUtil_snprintf( sTimeCoursePathAndStem, 
			  sizeof(sTimeCoursePathAndStem),
			  "%s/%s",
			  argv[nCurrentArg], argv[nCurrentArg+1] );
	  bLoadingTimeCourse = TRUE;
	  nCurrentArg += 2;
	}
  
      } else if( MATCH( sArg, "-overlay" ) ) {
	
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy arg into a destructible string */
	  DebugNote( ("Parsing -overlay option") );
	  xUtil_strncpy( sOverlayPathAndStem, argv[nCurrentArg+1],
			 sizeof(sOverlayPathAndStem) );
	  bLoadingOverlay = TRUE;
	  nCurrentArg += 2;
	  
	} else {
	  
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -overlay option",
			    "Expected an argument",
			    "This option needs an argument, the path and "
			    "stem of the data." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-overlay-reg" ) ) {
  
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy arg  */
	  DebugNote( ("Parsing -overlay-reg option") );
	  psOverlayRegistration = sOverlayRegistration;
	  xUtil_strncpy( sOverlayRegistration, argv[nCurrentArg+1],
			 sizeof(sOverlayRegistration) );
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -overlay-reg option",
			    "Expected an argument",
			    "This option needs an argument, the file name "
			    "of the registration data." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-overlay-offset" ) ) {
  
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy arg  */
	  DebugNote( ("Parsing -overlay-offset option") );
	  psOverlayOffsetPathAndStem = sOverlayOffsetPathAndStem;
	  xUtil_strncpy( sOverlayOffsetPathAndStem, argv[nCurrentArg+1],
			 sizeof(sOverlayOffsetPathAndStem) );
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -overlay-offset option",
			    "Expected an argument",
			    "This option needs an argument, the path "
			    "and stem of the offset volume." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-timecourse" ) ) {
  
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 ) {
	  
	  /* copy arg into a destructible string */
	  DebugNote( ("Parsing -timecourse option") );
	  xUtil_strncpy( sTimeCoursePathAndStem, argv[nCurrentArg+1],
			 sizeof(sTimeCoursePathAndStem) );
	  bLoadingTimeCourse = TRUE;
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -timecourse option",
			    "Expected an argument",
			    "This option needs an argument, the path and "
			    "stem of the data." );
	  nCurrentArg ++;
	}
	
      } else if( MATCH( sArg, "-timecourse-reg" ) ) {
  
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy arg  */
	  DebugNote( ("Parsing -timecourse-reg option") );
	  psTimeCourseRegistration = sTimeCourseRegistration;
	  xUtil_strncpy( sTimeCourseRegistration, argv[nCurrentArg+1],
			 sizeof(sTimeCourseRegistration) );
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -timecourse-reg option",
			    "Expected an argument",
			    "This option needs an argument, the file name "
			    "of the registration data." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-timecourse-offset" ) ) {
  
	/* make sure there are enough args */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy arg  */
	  DebugNote( ("Parsing -timecourse-offset option") );
	  psTimeCourseOffsetPathAndStem = sTimeCourseOffsetPathAndStem;
	  xUtil_strncpy( sTimeCourseOffsetPathAndStem, argv[nCurrentArg+1],
			 sizeof(sTimeCourseOffsetPathAndStem) );
	  nCurrentArg += 2;
	  
	} else {
    
	  /* misuse of that option */
	  tkm_DisplayError( "Parsing -timecourse-offset option",
			    "Expected an argument",
			    "This option needs an argument, the path "
			    "and stem to the offset volume." );
	  nCurrentArg ++;
	}
  
      } else if( MATCH( sArg, "-register" ) ) {
	
	/* set our flag */
	DebugNote( ("Enabling registration.") );
	bEnablingRegistration = TRUE;
	nCurrentArg ++;
	
      } else if( MATCH( sArg, "-segmentation" ) ||
		 MATCH( sArg, "-seg" ) ||
		 MATCH( sArg, "-parcellation" ) ||
		 MATCH( sArg, "-parc" ) ) {
	
	/* make sure there are enough args */
	if( argc > nCurrentArg + 2 &&
	    '-' != argv[nCurrentArg+1][0] &&
	    '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* copy path and color file */
	  DebugNote( ("Parsing -segmentation option") );
	  xUtil_strncpy( sSegmentationPath, argv[nCurrentArg+1],
			 sizeof(sSegmentationPath) );
	  xUtil_strncpy( sSegmentationColorFile, argv[nCurrentArg+2],
			 sizeof(sSegmentationColorFile) );
	  bLoadingSegmentation = TRUE;
	  nCurrentArg += 3;
	  
	} else {
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -segmentation/seg option",
			    "Expected two arguments",
			    "This option needs two arguments: the path of "
			    "the COR volume and the name of the colors "
			    "file." );
	  nCurrentArg ++;
	}
	
      } else if( MATCH( sArg, "-aux-segmentation" ) ||
		 MATCH( sArg, "-aux-seg" ) ) {
	
	/* make sure there are enough args */
	if( argc > nCurrentArg + 2 &&
	    '-' != argv[nCurrentArg+1][0] &&
	    '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* copy path and color file */
	  DebugNote( ("Parsing -aux-segmentation option") );
	  xUtil_strncpy( sAuxSegmentationPath, argv[nCurrentArg+1],
			 sizeof(sAuxSegmentationPath) );
	  xUtil_strncpy( sAuxSegmentationColorFile, argv[nCurrentArg+2],
			 sizeof(sAuxSegmentationColorFile) );
	  bLoadingAuxSegmentation = TRUE;
	  nCurrentArg += 3;
	  
	} else {
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -segmentation/seg option",
			    "Expected two arguments",
			    "This option needs two arguments: the path of "
			    "the COR volume and the name of the colors "
			    "file." );
	  nCurrentArg ++;
	}
	
      } else if( MATCH( sArg, "-voxel-label" ) ) {
	
	/* make sure there are enough args */
	if( argc > nCurrentArg + 2 &&
	    '-' != argv[nCurrentArg+1][0] &&
	    '-' != argv[nCurrentArg+2][0] ) {
	  
	  /* copy both file names */
	  DebugNote( ("Parsing -voxel-label option") );
	  xUtil_strncpy( sVLIFile1, argv[nCurrentArg+1], sizeof(sVLIFile1) );
	  xUtil_strncpy( sVLIFile2, argv[nCurrentArg+2], sizeof(sVLIFile2) );
	  bLoadingVLI = TRUE;
	  nCurrentArg += 3;
	  
	} else {
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -voxel-label option",
			    "Expected two arguments",
			    "This option needs two arguments: the file "
			    "names of the VLI volumes." );
	  nCurrentArg ++;
	}
	
      } else if( MATCH( sArg, "-f" ) ) {
  
	/* make sure subject is not already declared */
	if( bSubjectDeclared ) {
	  tkm_DisplayError( "Parsing -f option",
			    "Subject already declared",
			    "The -f option is only to be used if the "
			    "subject is not to be declared using the "
			    "subect/image type format." );
	  bFatalError = TRUE;
	  
	  /* check for path */
	} else if( argc > nCurrentArg + 1 &&
		   '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* read the path */
	  DebugNote( ("Parsing -f option") );
	  xUtil_strncpy( sSubject, argv[nCurrentArg+1], sizeof(sSubject) );
	  bUsingMRIRead = TRUE;
	  bSubjectDeclared = TRUE;
	  nCurrentArg += 2;
    
	  /* save subject home */
	  DebugNote( ("Setting subject home directory to %s", sSubject) );
	  SetSubjectHomeDir( sSubject );
	  
	  /* disable automatic file name making because there's no
	     home dir, really */
	  gEnableFileNameGuessing = FALSE;
	  
	} else {
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -f option",
			    "Expected an argument",
			    "This option needs an argument: the path or "
			    "file name of the data to use." );
	  nCurrentArg ++;
	}
	
      } else if( MATCH( sArg, "-fthresh" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -fthresh option") );
	  min = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
	  bThresh = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -fthresh option",
			    "Expected an argument",
			    "This option needs an argument: the threshold "
			    "value to use." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-fmid" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -fmid option") );
	  mid = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
	  bMid = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -fmid option",
			    "Expected an argument",
			    "This option needs an argument: the threshold "
			    "value to use." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-sdir" ) ) {
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -sdir option") );
	  gsCommandLineSubjectsDir = argv[nCurrentArg+1] ;
	  nCurrentArg +=2 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -fslope option",
			    "Expected an argument",
			    "This option needs an argument: the threshold "
			    "value to use." );
	  nCurrentArg += 1;
	}
      } else if( MATCH( sArg, "-fslope" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -fslope option") );
	  slope = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
	  bSlope = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -fslope option",
			    "Expected an argument",
			    "This option needs an argument: the threshold "
			    "value to use." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-fsmooth" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -fsmooth option") );
	  smoothSigma = (float) atof( argv[nCurrentArg+1] );
	  bSmooth = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -fsmooth option",
			    "Expected an argument",
			    "This option needs an argument: the sigma of the Gaussian "
			    "value to use." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-revphaseflag" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -revphaseflag option") );
	  nRevPhaseFlag = atoi( argv[nCurrentArg+1] );
	  bRevPhaseFlag = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -revphaseflag option",
			    "Expected an argument",
			    "This option needs an argument: a 1 or 0 to "
			    "turn it on or off." );
	  nCurrentArg += 1;
	}
	
      } else if( MATCH( sArg, "-truncphaseflag" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -truncphaseflag option") );
	  nTruncPhaseFlag = atoi( argv[nCurrentArg+1] );
	  bTruncPhaseFlag = TRUE;
	  nCurrentArg +=2 ;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -truncphaseflag option",
			    "Expected an argument",
			    "This option needs an argument: a 1 or 0 to "
			    "turn it on or off." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-segmentation-opacity" ) ||
		 MATCH( sArg, "-roialpha" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1 &&
	    '-' != argv[nCurrentArg+1][0] ) {
    
	  /* get the value */
	  DebugNote( ("Parsing -roialpha option") );
	  fSegmentationAlpha = (float) atof( argv[nCurrentArg+1] );
	  bSegmentationAlpha = TRUE;
	  nCurrentArg +=2 ;
    
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -roialpha option",
			    "Expected an argument",
			    "This option needs an argument: the alpha "
			    "value to use." );
	  nCurrentArg += 1;
	}
	
      } else if( MATCH( sArg, "-aux" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* read in the aux file name */
	  DebugNote( ("Parsing -aux option") );
	  xUtil_strncpy( sAuxVolume, argv[nCurrentArg+1], 
			 sizeof(sAuxVolume) );
	  bLoadingAuxVolume = TRUE;
	  nCurrentArg += 2;
	  
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -aux option",
			    "Expected an argument",
			    "This option needs an argument: the image type, "
			    "directory, or file name of the data to load "
			    "as the aux volume." );
	  nCurrentArg += 1;
	}     
  
      } else if( MATCH( sArg, "-main-transform" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* read in the main transform file name */
	  DebugNote( ("Parsing -main-transform option") );
	  xUtil_strncpy( sMainTransform, argv[nCurrentArg+1], 
			 sizeof(sMainTransform) );
	  bLoadingMainTransform = TRUE;
	  nCurrentArg += 2;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -main-transform option",
			    "Expected an argument",
			    "This option needs an argument: the "
			    "file name of the transform to load." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-aux-transform" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* read in the aux transform file name */
	  DebugNote( ("Parsing -aux-transform option") );
	  xUtil_strncpy( sAuxTransform, argv[nCurrentArg+1], 
			 sizeof(sAuxTransform) );
	  bLoadingAuxTransform = TRUE;
	  nCurrentArg += 2;
	  
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -aux-transform option",
			    "Expected an argument",
			    "This option needs an argument: the "
			    "file name of the transform to load." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-label" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
    
	  /* read in the label name */
	  DebugNote( ("Parsing -label option") );
	  xUtil_strncpy( sLabel, argv[nCurrentArg+1], sizeof(sLabel) );
	  
	  bLoadingLabel = TRUE;
	  nCurrentArg += 2;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -label option",
			    "Expected an argument",
			    "This option needs an argument: the file name "
			    "of the label file." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-headpts" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* read in the head pts and transform file name */
	  DebugNote( ("Parsing -headpts option") );
	  xUtil_strncpy( sHeadPts, argv[nCurrentArg+1], sizeof(sHeadPts) );
	  
	  /* if they gave us a transform file as well... */
	  if( argc > nCurrentArg + 2
	      && '-' != argv[nCurrentArg+2][0] ) {
	    
	    /* save that */
	    DebugNote( ("Parsing transform of -headpts option") );
	    xUtil_strncpy( sHeadPtsTransform, argv[nCurrentArg+2],
			   sizeof(sHeadPtsTransform));
	    bHaveHeadPtsTransform = TRUE;
	  } else {
	    bHaveHeadPtsTransform = FALSE;
	  }
	  
	  bLoadingHeadPts = TRUE;
	  nCurrentArg += 3;
	  
	} else { 
	  
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -headpts option",
			    "Expected one or two arguments",
			    "This option needs an argument: the file name "
			    "of the head points file, and optionally the "
			    "file name of the transform to use." );
	  nCurrentArg += 1;
	}
	
      } else if( MATCH( sArg, "-overlaycache" ) ) {
	
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* get the value */
	  DebugNote( ("Parsing -overlaycache option") );
	  nUseOverlayCacheFlag = atoi( argv[nCurrentArg+1] );
	  bUseOverlayCacheFlag = TRUE;
	  nCurrentArg +=2 ;
    
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -overlaycache option",
			    "Expected an argument",
			    "This option needs an argument: a 1 or 0 to "
			    "turn the option on or off." );
	  nCurrentArg += 1;
	}
  
	/* rkt - commented out because the functional volume should no
	   longer set the conversion method explicitly. it should only be
	   set when parsing the register.dat file. */
#if 0
      } else if( MATCH( sArg, "-float2int" ) ) {
  
	/* check for the value following the switch */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
    
	  /* get the value */
	  DebugNote( ("Parsing -float2int option") );
	  if( MATCH( argv[nCurrentArg+1], "tkreg" ) ) {
	    convMethod = FunD_tConversionMethod_FCF;
	    bSetConversionMethod = TRUE;
	    nCurrentArg +=2;
	  }
	  else if( MATCH( argv[nCurrentArg+1], "floor" ) ) {
	    convMethod = FunD_tConversionMethod_FFF;
	    bSetConversionMethod = TRUE;
	    nCurrentArg +=2;
	  }
	  else if( MATCH( argv[nCurrentArg+1], "round" ) ) {
	    convMethod = FunD_tConversionMethod_Round;
	    bSetConversionMethod = TRUE;
	    nCurrentArg +=2;
	  } else {
	    tkm_DisplayError( "Parsing -float2int option",
			      "Argument not recognized",
			      "Please specify tkreg, floor, or round "
			      "as the conversion method." );
	    nCurrentArg +=1;
	  }
    
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -overlaycache option",
			    "Expected an argument",
			    "This option needs an argument: a 1 or 0 to "
			    "turn the option on or off." );
	  nCurrentArg += 1;
	}
#endif
  
      } else if( MATCH( sArg, "-interface" ) ) {
  
	/* check for another value */
	if( argc > nCurrentArg + 1
	    && '-' != argv[nCurrentArg+1][0] ) {
	  
	  /* copy in the interface name */
	  DebugNote( ("Parsing -interface option") );
	  strcpy( gInterfaceScriptName, argv[nCurrentArg+1] );
	  nCurrentArg += 2;
	  
	} else { 
    
	  /* misuse of that switch */
	  tkm_DisplayError( "Parsing -interface option",
			    "Expected an argument",
			    "This option needs an argument: the file name "
			    "of the interface script file to use." );
	  nCurrentArg += 1;
	}
  
      } else if( MATCH( sArg, "-" ) ) {
  
	/* set no edit mode. */
	DebugNote( ("Parsing - option") );
	bNoEdit = TRUE;
	nCurrentArg += 1;
	
      } else if( MATCH( sArg, "-csurf" ) ) {
	
	/* set name of interface to tkmedit_csurf.tcl */
	DebugNote( ("Parsing -csurf option") );
	gbUseCsurfInterface = TRUE;
	nCurrentArg ++;
  
      } else {
  
	/* unrecognized option, build an error message and ignore it. */
	xUtil_snprintf( sError, sizeof(sError), 
			"Option %s not recognized", sArg );
	tkm_DisplayError( "Parsing command line options",
			  sError,
			  "This option was not recognized and ignored." );
	nCurrentArg ++;
      }
      
      /* check for local keyword */
    } else if ( MATCH( sArg, "local" ) || 
		MATCH( sArg, "." ) ) {
      
      /* make sure subject is not already declared */
      if( bSubjectDeclared ) {
	tkm_DisplayError( "Parsing local option",
			  "Subject already declared",
			  "The local option is only to be used if the "
			  "subject is not to be declared using the "
			  "subect/image type format." );
	bFatalError = TRUE;
  
      } else {
  
	/* set local flag */
	DebugNote( ("Parsing local option") );
	bLocalImageDir = TRUE;
	bSubjectDeclared = TRUE;
	nCurrentArg ++;
	
	/* save subject home */
	DebugNote( ("Setting user home dir to %s", gsUserHomeDir) );
	SetSubjectHomeDir( gsUserHomeDir );
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
	    tkm_DisplayError( "Parsing subject and image name",
			      "Image name not declared",
			      "When declaring the subject and image type, "
			      "two arguments are expected." );
	    bFatalError = TRUE;
	    
	  } else { 
	    
	    /* read in subject and image name. */
	    DebugNote( ("Parsing subject name") );
	    xUtil_strncpy( sSubject, argv[nCurrentArg], sizeof(sSubject) );
	    
	    DebugNote( ("Parsing image type") );
	    xUtil_strncpy( sImageDir, argv[nCurrentArg+1], 
			   sizeof(sImageDir) );
	    bSubjectDeclared = TRUE;
	    nCurrentArg += 2;
	    
	    /* save subject home */
	    DebugNote( ("Setting subject home from env") );
	    eResult = SetSubjectHomeDirFromEnv( sSubject );
	    DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
	  }
	  
	  /* check for a surface. if we have enough args... */
	  if( argc > nCurrentArg ) {
	    
	    /* and this one isn't a switch or the local flag... */
	    if( '-' != argv[nCurrentArg][0] 
		&& !MATCH( argv[nCurrentArg], "local" )
		&& !MATCH( argv[nCurrentArg], "." ) ) {
	      
	      /* read it as a surface. */
	      DebugNote( ("Parsing surface name") );
	      xUtil_strncpy( sSurface, argv[nCurrentArg], sizeof(sSurface) );
	      bSurfaceDeclared = TRUE;
	      nCurrentArg ++;
	    }
	  }
	} else {
	  tkm_DisplayError( "Parsing arguments",
			    "Subject and image not declared",
			    "You must pass a subject name and image type "
			    "to load. Enter tkmedit with no arguments "
			    "for a list of parameters." );
	  bFatalError = TRUE;
	}
  
      } else {
  
	/* totally unrecognized */
	xUtil_snprintf( sError, sizeof(sError),
			"Unrecognized option: %s", sArg );
	tkm_DisplayError( "Parsing arguments",
			  sError,
			  "This option is not recognized by tkmedit." );
	nCurrentArg++;
      }
    }
  }
  
  /* check for fatal error. */
  if( bFatalError ) {
    /* probably have some error messages that didn't get printed cuz
       tcl wasn't loaded yet, so flush them to shell now. */
    PrintCachedTclErrorDlogsToShell();
    DebugPrint( ( "Fatal error in parsing command line args.\n" ) );
    exit(1);
  }
  
  /* if using local directory, copy 'local' in subject name for
     historical reasons. */
  if( bLocalImageDir ) {
    DebugNote( ("Copying local into subject name") );
    xUtil_strncpy( sSubject, "local", sizeof(sSubject) );
    DebugNote( ("Copying empty string into image dir") );
    xUtil_strncpy( sImageDir, "", sizeof(sImageDir) );
  }
  
  /* disable editing */
  if( bNoEdit ) {
    editflag = FALSE;
  }
  
  /* bleah. if we have the using mri read flag set, load the images from
     the subject name, otherwise use the image dir. */
  if( bUsingMRIRead ) {
    DebugNote( ("Loading volume %s", sSubject) );
    eResult = LoadVolume( tkm_tVolumeType_Main, sSubject );
    if( tkm_tErr_NoErr != eResult ) {
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }
  } else {
    DebugNote( ("Loading volume %s", sImageDir) );
    eResult = LoadVolume( tkm_tVolumeType_Main, sImageDir );
    if( tkm_tErr_NoErr != eResult ) {
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }
    
    /* check to see if we don't have a subject */
    Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main],
			  sSubjectTest, sizeof(sSubjectTest) );
    if( strcmp( sSubjectTest, "" ) == 0 ) {
      /* manually set the subject and image name */
      Volm_SetSubjectName( gAnatomicalVolume[tkm_tVolumeType_Main], 
			   sSubject );
      Volm_SetVolumeName( gAnatomicalVolume[tkm_tVolumeType_Main], 
			  sImageDir );
    }
  }
  
  /* If we got a non-default brightness and contrast, set it now. */
  if( bBrightContrastMain ) {
    SetVolumeBrightnessAndContrast( tkm_tVolumeType_Main, 
				    fBrightnessMain, fContrastMain );
  }

  /* if reading in an aux image... */
  if( bLoadingAuxVolume ) {
    DebugNote( ("Loading aux volume %s", sAuxVolume) );
    eResult = LoadVolume( tkm_tVolumeType_Aux, sAuxVolume );

    /* If we got a non-default brightness and contrast, set it now. */
    if( bBrightContrastAux ) {
      SetVolumeBrightnessAndContrast( tkm_tVolumeType_Aux, 
				      fBrightnessAux, fContrastAux );
    }
  }
  
  /* load in the display transforms. */
  if( bLoadingMainTransform ) {
    DebugNote( ("Loading main display transform %s", sMainTransform) );
    eResult = LoadDisplayTransform( tkm_tVolumeType_Main, sMainTransform );
  }
  if( bLoadingAuxTransform ) {
    DebugNote( ("Loading aux display transform %s", sAuxTransform) );
    eResult = LoadDisplayTransform( tkm_tVolumeType_Aux, sAuxTransform );
  }
  
  /* load surface. trasnsform must be inited first. */
  if ( bSurfaceDeclared ) {
    DebugNote( ("Loading surface") );
    eResult = LoadSurface( tkm_tSurfaceType_Main, sSurface );
  }
  
  /* load segmentation */
  if( bLoadingSegmentation ) {
    eResult = LoadSegmentationVolume( tkm_tSegType_Main, sSegmentationPath,
				      sSegmentationColorFile );
    /* set roi alpha */
    if( bSegmentationAlpha ) {
      SetSegmentationAlpha( fSegmentationAlpha );
    }
  }
  
  /* load aux segmentation */
  if( bLoadingAuxSegmentation ) {
    eResult = LoadSegmentationVolume( tkm_tSegType_Aux, sAuxSegmentationPath,
				      sAuxSegmentationColorFile );
  }
  
  /* load VLIs */
  if( bLoadingVLI ) {
    eResult = LoadVLIs( sVLIFile1, sVLIFile2 );
  }

  /* load the label */
  if( bLoadingLabel ) {
    eResult = LoadSelectionFromLabelFile( sLabel );
  }
  
  /* load head pts */
  if( bLoadingHeadPts ) {
    if( bHaveHeadPtsTransform ) {
      eResult = LoadHeadPts( sHeadPts, sHeadPtsTransform ); 
    } else {
      eResult = LoadHeadPts( sHeadPts, NULL );
    }
  }
  
  /* load functional overlay data */
  if( bLoadingOverlay ) {
    eResult = LoadFunctionalOverlay( sOverlayPathAndStem,
             psOverlayOffsetPathAndStem,
             psOverlayRegistration );
    
    if( eResult == tkm_tErr_NoErr && bSmooth ) {
      eResult = SmoothOverlayData( smoothSigma );
    }
  }
  
  /* load functional time course data */
  if( bLoadingTimeCourse ) {
    eResult = LoadFunctionalTimeCourse( sTimeCoursePathAndStem, 
          psTimeCourseOffsetPathAndStem,
          psTimeCourseRegistration );
  }
  
  /* set registration */
  DebugNote( ("%sabling registration", bEnablingRegistration?"En":"Dis") );
  FunV_EnableRegistration( gFunctionalVolume, bEnablingRegistration );
  
  /* if regisration is enabled, set conversion method to round. */
  if( bEnablingRegistration ) {
    FunV_SetConversionMethod( gFunctionalVolume, 
            FunD_tConversionMethod_Round );
  }
  
  
  /* set functional color scale stuff */
  if( bThresh || bMid || bThresh ) {
    DebugNote( ("Setting functional threshold min=%f mid=%f slope=%f\n",
    min, mid, slope ) );
    eFunctional = FunV_SetThreshold( gFunctionalVolume, min, mid, slope );
  }
  if( bTruncPhaseFlag ) {
    DebugNote( ("Setting trunc phase flag to %d\n", nTruncPhaseFlag) );
    eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
               FunV_tDisplayFlag_Ol_TruncateNegative,
               (tBoolean) nTruncPhaseFlag );
  }
  if( bRevPhaseFlag ) {
    DebugNote( ("Setting rev phase flag to %d\n", nRevPhaseFlag ) );
    eFunctional = FunV_SetDisplayFlag( gFunctionalVolume, 
               FunV_tDisplayFlag_Ol_ReversePhase,
               (tBoolean) nRevPhaseFlag );
  }
  if( bUseOverlayCacheFlag ) { 
    DebugNote( ("Setting overlay cache flag to %d\n", 
    nUseOverlayCacheFlag) );
    eFunctional = FunV_UseOverlayCache( gFunctionalVolume,
          (tBoolean) nUseOverlayCacheFlag );
  }
  
  /* rkt - commented out because the functional volume should no
     longer set the conversion method explicitly. it should only be
     set when parsing the register.dat file. */
#if 0
  if( bSetConversionMethod ) {
    DebugNote( ("Setting conversion method to %d", convMethod) );
    eFunctional = FunV_SetConversionMethod( gFunctionalVolume, convMethod );
  }
#endif
  
  /* clear error flag because if we get here we've already handled it
     with an error message. */
  eResult = tkm_tErr_NoErr;
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void save_rgb( char* fname ) {
  
  GLint rowlength, skiprows, skippixels, alignment;
  GLboolean swapbytes, lsbfirst;
  int nNumBytes = 0;
  GLenum eGL;
  RGB_IMAGE *image;
  int y,xloc,yloc,width,height,size;
  unsigned short *r,*g,*b;
  unsigned short  *red = NULL;
  unsigned short  *green = NULL;
  unsigned short  *blue = NULL;
  FILE *fp;
  
  if( NULL == gMeditWindow )
    return;
  
  MWin_GetWindowSize( gMeditWindow,
          &xloc, &yloc, &width, &height );
  
  size = width*height;
  nNumBytes = sizeof( unsigned short ) * size;
  
  red  = (unsigned short*) malloc( nNumBytes );
  if( NULL == red ) {
    DebugPrint( ( "save_rgb: allocation of red buffer failed.\n" ) );
    goto cleanup;
  }
  green = (unsigned short*) malloc( nNumBytes );
  if( NULL == green ) {
    DebugPrint( ( "save_rgb: allocation of green buffer failed.\n" ) );
    goto cleanup;
  }
  blue  = (unsigned short*) malloc( nNumBytes );
  if( NULL == blue ) {
    DebugPrint( ( "save_rgb: allocation of blue buffer failed.\n" ) );
    goto cleanup;
  }
  
  glReadBuffer( GL_FRONT );
  
  glGetBooleanv(GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetBooleanv(GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment);
  
  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
  
  glReadPixels(0, 0, width, height, GL_RED,  
         GL_UNSIGNED_SHORT, (GLvoid *)red);
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL ) {
    DebugPrint( ( "glReadPixels got error %d\n", eGL ) );
    goto cleanup;
  }
  glReadPixels(0, 0, width, height, GL_GREEN,
         GL_UNSIGNED_SHORT, (GLvoid *)green);
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL ) {
    DebugPrint( ( "glReadPixels got error %d\n", eGL ) );
    goto cleanup;
  }
  glReadPixels(0, 0, width, height, GL_BLUE, 
         GL_UNSIGNED_SHORT, (GLvoid *)blue);
  eGL = glGetError ();
  if( GL_NO_ERROR != eGL ) {
    DebugPrint( ( "glReadPixels got error %d\n", eGL ) );
    goto cleanup;
  }
  
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
  
 cleanup:
  
  if( NULL != red ) 
    free( red ); 
  if( NULL != green ) 
    free( green ); 
  if( NULL != blue ) 
    free( blue );
}


void GotoSurfaceVertex ( Surf_tVertexSet iSurface, int inVertex ) {
  
  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    anaIdx;
  char      sDescription[STRLEN];
  char      sSetName[STRLEN];
  
  /* make sure we have a surface. */
  if( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the vertex */
  eSurface = Surf_GetNthVertex( gSurface[tkm_tSurfaceType_Main],
				iSurface, inVertex, &anaIdx,
				sDescription );
  if( Surf_tErr_NoErr != eSurface ) 
    goto error;
  
  /* print the result string */
  Surf_GetSurfaceSetName( iSurface, sSetName );
  OutputPrint "%s vertex index %d:\n\t%s\n", 
    sSetName, inVertex, sDescription EndOutputPrint;
  
  /* RKT: We don't need to adjust surface verts any more. */
#if 0
  /* adjust it so it aligns to the surface. */
  eWindow = MWin_AdjustSurfaceAnaIdx( gMeditWindow, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
#endif
  
  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  /* zoom around new cursor. */
  eWindow = MWin_SetZoomCenterToCursor ( gMeditWindow, -1 );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;
  
 error:
  
  DebugPrint( ( "Error in GotoSurfaceVertex( %d, %d )\n",
    (int)iSurface, inVertex ) );
  
 cleanup:
  return;
}

void FindNearestSurfaceVertex ( Surf_tVertexSet iSet ) {
  
  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    cursor;
  xVoxel    anaIdx;
  char      sDescription[STRLEN];
  char      sSetName[STRLEN];
  
  /* make sure we have a surface. */
  if( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the cursor */
  eWindow = MWin_GetCursor ( gMeditWindow, &cursor );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  /* RKT: We don't need to adjust surface verts any more. */
#if 0
  /* first unadjust the point. */
  eWindow = MWin_UnadjustSurfaceAnaIdx( gMeditWindow, &cursor );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
#endif
  
  /* get the vertex */
  eSurface = Surf_GetClosestVertexVoxel( gSurface[tkm_tSurfaceType_Main],
					 iSet, &cursor, &anaIdx,
					 sDescription);
  if( Surf_tErr_NoErr != eSurface ) 
    goto error;
  
  /* print the result string */
  Surf_GetSurfaceSetName( iSet, sSetName );
  OutputPrint "Nearest %s vertex to %d, %d, %d:\n\t%s\n",
    sSetName, xVoxl_ExpandInt( &cursor ), sDescription EndOutputPrint;
  
  /* RKT: We don't need to adjust surface verts any more. */
#if 0
  /* adjust it so it aligns to the surface. */
  eWindow = MWin_AdjustSurfaceAnaIdx( gMeditWindow, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
#endif
  
  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  /* zoom around new cursor. */
  eWindow = MWin_SetZoomCenterToCursor ( gMeditWindow, -1 );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;
  
 error:
  
  DebugPrint( ( "Error in FindNearestSurfaceVertex( %d )\n",
    (int)iSet ) );
  
 cleanup:
  return;
}

void FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet iSet ) {
  
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    cursor;
  xVoxel    origAnaIdx;
  xVoxel    interpAnaIdx;
  char      sDescription[STRLEN];
  char      sSetName[STRLEN];
  
  /* make sure we have a surface. */
  if( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the cursor */
  eWindow = MWin_GetCursor ( gMeditWindow, &cursor );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  /* RKT: We don't need to adjust surface verts any more. */
#if 0
  /* first unadjust the point. */
  eWindow = MWin_UnadjustSurfaceAnaIdx( gMeditWindow, &cursor );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
#endif
  
  /* get the verteices */
  eWindow = MWin_GetClosestInterpSurfVoxel( gMeditWindow,
					    tkm_tSurfaceType_Main,
					    iSet, &cursor, 
					    &origAnaIdx, &interpAnaIdx,
					    sDescription);
  if( MWin_tErr_NoErr != eWindow ) 
    goto error;
  
  /* print the result string */
  Surf_GetSurfaceSetName( iSet, sSetName );
  OutputPrint "Nearest %s vertex to %d, %d, %d:\n\t%s\n",
    sSetName, xVoxl_ExpandInt( &cursor ), sDescription EndOutputPrint;
  
  /* RKT: We don't need to adjust surface verts any more. */
#if 0
  /* adjust it so it aligns to the surface. */
  eWindow = MWin_AdjustSurfaceAnaIdx( gMeditWindow, &interpAnaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
#endif
  
  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &interpAnaIdx );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  /* zoom around new cursor. */
  eWindow = MWin_SetZoomCenterToCursor ( gMeditWindow, -1 );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;
  
 error:
  
  DebugPrint( ( "Error in FindNearestInterpolatedSurfaceVertex( %d )\n",
		(int)iSet ) );
  
 cleanup:
  return;
}

void AverageSurfaceVertexPositions ( int inNumAverages ) {
  
  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;

  /* make sure we have a surface. */
  if( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* call the average function. */
  eSurface = Surf_AverageVertexPositions( gSurface[tkm_tSurfaceType_Main],
					  inNumAverages );
  if( Surf_tErr_NoErr != eSurface )
    goto error;
  
  /* redraw the window. */
  eWindow = MWin_RedrawAll( gMeditWindow );
  if( MWin_tErr_NoErr != eWindow )
    goto error;
  
  goto cleanup;
  
 error:
  
  DebugPrint( ( "Error in AverageSurfaceVertexPositions( %d )\n",
		(int)inNumAverages ) );
  
 cleanup:
  return;
}

void WriteVoxelToControlFile ( xVoxelRef iAnaIdx ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen] = "";
  FILE*      file         = NULL;
  xVoxel    ras;
  
  DebugEnterFunction( ("WriteVoxelToControlFile ( iAnaIdx=%d,%d,%d )",
           xVoxl_ExpandInt( iAnaIdx )) );
  
  /* make the file name */
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, 
    sFileName, sizeof(sFileName) );
  
  /* open it */
  DebugNote( ("Opening control point file") );
  file = fopen( sFileName, "a+" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );
  
  /* convert idx to ras */
  eVolume = Volm_ConvertIdxToRAS( gAnatomicalVolume[tkm_tVolumeType_Main],
          iAnaIdx, &ras );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* write RAS space pt to file */
  DebugNote( ("Writing control point %.2f,%.2f,%.2f to file",
        xVoxl_ExpandFloat( &ras ) ));
  fprintf( file,"%f %f %f\n", xVoxl_ExpandFloat( &ras ) );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != file ) {
    DebugNote( ("Closing control point file") );
    fclose( file );
  }
  
  DebugExitFunction;
}

void WriteVoxelToEditFile ( xVoxelRef iAnaIdx ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen] = "";
  FILE*      file         = NULL;
  xVoxel    ras;
  xVoxel    tal;
  
  DebugEnterFunction( ("WriteVoxelToEditFile ( iAnaIdx=%d,%d,%d )",
           xVoxl_ExpandInt( iAnaIdx )) );
  
  /* make the file name */
  MakeFileName( "edit.dat", tkm_tFileName_Edit, 
    sFileName, sizeof(sFileName) );
  
  /* open it */
  DebugNote( ("Opening edit file") );
  file = fopen( sFileName, "w" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );
  
  /* convert idx to ras */
  eVolume = Volm_ConvertIdxToRAS( gAnatomicalVolume[tkm_tVolumeType_Main],
          iAnaIdx, &ras );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* write RAS space pt to file */
  DebugNote( ("Writing ras edit point %.2f,%.2f,%.2f to file",
        xVoxl_ExpandFloat( &ras ) ));
  fprintf( file,"%f %f %f\n", xVoxl_ExpandFloat( &ras ) );
  
  /* convert to tal and write that. */
  eVolume = Volm_ConvertIdxToTal( gAnatomicalVolume[tkm_tVolumeType_Main],
          iAnaIdx, &tal );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  DebugNote( ("Writing tal edit point %.2f,%.2f,%.2f to file",
        xVoxl_ExpandFloat( &tal ) ));
  fprintf( file,"%f %f %f\n", xVoxl_ExpandFloat( &tal ) );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != file ) {
    DebugNote( ("Closing edit file") );
    fclose( file );
  }
  
  DebugExitFunction;
}

void GotoSavedCursor () {
  
  char sFileName[tkm_knPathLen] = "";
  
  DebugEnterFunction( ("GotoSavedCursor ()") );
  
  /* make the file name. */
  DebugNote( ("Making file name from edit.dat") );
  MakeFileName( "edit.dat", tkm_tFileName_Edit, sFileName, sizeof(sFileName) );
  
  ReadCursorFromEditFile( sFileName );
  
  DebugExitFunction;
}

void GotoAnotherSubjectSavedCursor ( char* isSubjectName ) {
  
}

void ReadCursorFromEditFile ( char* isFileName ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  FILE*      file    = NULL;
  float      fRASX   = 0;
  float      fRASY   = 0;
  float      fRASZ   = 0;
  xVoxel    ras;
  xVoxel    idx;
  
  DebugEnterFunction( ("ReadCursorFromEditFile ( isFileName=%s )",
           isFileName) );
  
  /* open it. */
  DebugNote( ("Opening edit file") );
  file = fopen( isFileName, "r" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );
  
  /* read the point */
  DebugNote( ("Reading point from edit file") );
  fscanf( file, "%f %f %f", &fRASX, &fRASY, &fRASZ );
  xVoxl_SetFloat( &ras, fRASX, fRASY, fRASZ );
  
  /* convert to volume voxel. */
  eVolume = Volm_ConvertRASToIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
          &ras, &idx );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* build and set cursor */
  DebugNote( ("Setting cursor in main window") );
  MWin_SetCursor( gMeditWindow, -1, &idx );
  DebugNote( ("Setting zoom center in main window") );
  MWin_SetZoomCenterToCursor( gMeditWindow, -1 );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  tkm_DisplayError( "Going to saved cursor",
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the cursor from the edit.dat "
        "file." );
  EndDebugCatch;
  
  if( NULL != file ) {
    DebugNote( ("Closing edit file") );
    fclose( file );
  }
  
  DebugExitFunction;
}

void UpdateAndRedraw () {
  
  if ( NULL != gMeditWindow ) {
    MWin_RedrawAll( gMeditWindow );
  }
}


void tkm_HandleIdle () {
  
  /* just call the tk event handling function */
  Tk_DoOneEvent( TK_ALL_EVENTS | TK_DONT_WAIT );
}

// ================================================================== SURFACES

tkm_tErr LoadSurface ( tkm_tSurfaceType iType,
		       char*    isName ) {
  
  tkm_tErr  eResult           = tkm_tErr_NoErr;
  Surf_tErr eSurface           = Surf_tErr_NoErr;
  Trns_tErr eTrns             = Trns_tErr_NoErr;
  tBoolean  bLoaded           = FALSE;
  char      sName[tkm_knPathLen]       = "";
  char      sError[tkm_knErrStringLen] = "";
  mriTransformRef surfaceTransform = NULL;
#if 0
  MATRIX *tmp1, *tmp2;
#endif

  DebugEnterFunction( ("LoadSurface( iType=%d, isName=%s )", 
           (int)iType, isName) );
  
  DebugAssertThrowX( (NULL != isName),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumSurfaceTypes),
		     eResult, tkm_tErr_InvalidParameter );
  
  /* make file name */
  DebugNote( ("Making file name from %s", isName) );
  MakeFileName( isName, tkm_tFileName_Surface, sName, sizeof(sName) );
  
  /* delete existing if we already have it */
  if( NULL != gSurface[iType] ) {
    Surf_Delete( &gSurface[iType] );
  }


  DebugNote( ("Cloning gIdxToRASTransform to get surfaceTransform") );
  eTrns = Trns_DeepClone( gIdxToRASTransform, &surfaceTransform );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTrns),
		     eResult, tkm_tErr_CouldntAllocate );


  /* RKT: Why was this happening? This makes BtoRAS be _not_ the
     inverse of RAStoB, which is bad. It seems necessary to make it
     work, but this should be fixed some times soon. */

  // modify surfaceTransform->mBtoRAS
  *MATRIX_RELT(surfaceTransform->mBtoRAS, 1, 4) = 128;
  *MATRIX_RELT(surfaceTransform->mBtoRAS, 2, 4) = -128;
  *MATRIX_RELT(surfaceTransform->mBtoRAS, 3, 4) = 128;

  /* RKT: We don't need to do this any more because mriSurface only
     uses BtoRAS and BRAStoB. I don't know why this stuff was being
     done anyway. */

  // modify surfaceTransform->mARAStoBRAS
#if 0
  tmp1 = MatrixInverse(surfaceTransform->mAtoRAS, NULL);
  tmp2 = MatrixMultiply(surfaceTransform->mAtoB, tmp1, NULL);
  surfaceTransform->mARAStoBRAS = 
    MatrixMultiply(surfaceTransform->mBtoRAS, tmp2,
		   surfaceTransform->mARAStoBRAS);

#if 0
  DebugPrint(("AtoRAS-1\n"));
  MatrixPrint(stderr,tmp1);
  DebugPrint(("AtoB\n"));
  MatrixPrint(stderr,surfaceTransform->mAtoB);
  DebugPrint(("tmp2\n"));
  MatrixPrint(stderr,tmp2);
  DebugPrint(("BtoRAS\n"));
  MatrixPrint(stderr,surfaceTransform->mBtoRAS);
  DebugPrint(("ARAStoBRAS\n"));
  MatrixPrint(stderr,surfaceTransform->mARAStoBRAS);
#endif

  MatrixFree(&tmp1);
  MatrixFree(&tmp2);
#endif

  /* RKT: So at this point, AtoRAS is the same as it was in
     gIdxToRASTransform, which is extract_i_to_r. BtoRAS is almost the
     same as it was, except the transformation part of the matrix is
     128, -128, 128. ARAStoBRAS is all weird. BUT, much of this
     doesn't matter, since in mriSurface, only ConvertBtoRAS and
     ConvertBRAStoB is called on this transform, so we're only really
     using BtoRAS and BRAStoB. BRAStoB is still the calculated inverse
     of the original BtoRAS.

     Note that at this point, surfaceTransform is in an invalid
     state. Whenever Trns_CopyAtoRAS, Trns_CopyBtoRAS, or
     Trns_CopyARAStoBRAS are called, mriTransform automatically
     calculates the inverse and compositions of these matrices (namely
     ARAStoA, BRAStoB, BRAStoARAS, AtoB, and BtoA). But, since this
     code modified the member matrices directly, that automatic
     calculation has never been done. */


  /* create the surface */
  DebugNote( ("Creating surface") );
  eSurface = Surf_New( &gSurface[iType], sName, surfaceTransform);
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
         eResult, tkm_tErr_CouldntLoadSurface );


#if 0
  printf("LoadSurface surfaceTransform================================\n");
  Trns_DebugPrint_( surfaceTransform );
#endif
 
  /* see if it was loaded */
  DebugNote( ("Loading main vertex set") );
  eSurface = Surf_IsVertexSetLoaded( gSurface[iType], Surf_tVertexSet_Main, 
             &bLoaded );
  DebugAssertThrowX( (bLoaded), eResult, tkm_tErr_CouldntLoadSurface );
  
  /* set the medit window surface. */
  DebugNote( ("Setting surface in main window") );
  MWin_SetSurface( gMeditWindow, -1, iType, gSurface[iType] );
  
  /* turn on the loading and viewing options for surfaces in our interface.
     also turn the surface display onn in the window. turn on the 
     interpolated vertex display. */
  DebugNote( ("Setting tcl options") );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, 
		      bLoaded ? "1" : "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceViewingOptions, 
		      bLoaded ? "1" : "0" );
  MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_MainSurface, 
		       bLoaded ? TRUE : FALSE);
  MWin_SetDisplayFlag( gMeditWindow, -1,
		       DspA_tDisplayFlag_InterpolateSurfaceVertices,
		       bLoaded ? TRUE : FALSE);
  
  /* load other vertices. if these fail, it's okay. */
  DebugNote( ("Loading orig set") );
  LoadSurfaceVertexSet( iType, Surf_tVertexSet_Original, "orig" );
  DebugNote( ("Loading pial set") );
  LoadSurfaceVertexSet( iType, Surf_tVertexSet_Pial, "pial" );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading Surface %s", isName );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the surface you specified. "
        "This could be because the file wasn't a surface "
        "type that tkmedit recognizes "
        "or the file was unreadable due to permissions." );
  EndDebugCatch;
  
  Trns_Delete( &surfaceTransform );

  DebugExitFunction;
  
  return eResult;
}

tkm_tErr LoadSurfaceVertexSet ( tkm_tSurfaceType iType,
        Surf_tVertexSet   iSet,
        char*     isName ) {
  
  tkm_tErr      eResult  = tkm_tErr_NoErr;
  Surf_tErr      eSurface = Surf_tErr_NoErr;
  tBoolean      bLoaded  = FALSE;
  DspA_tDisplayFlag flag     = DspA_tDisplayFlag_None;
  tkm_tTclCommand   command  = 0;
  
  DebugEnterFunction( ("LoadSurfaceVertexSet( iType=%d, iSet=%d, isName=%s )", 
           (int)iType, (int)iSet, isName) );
  
  DebugAssertThrowX( (NULL != isName),
         eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumSurfaceTypes),
         eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != gSurface[iType]),
         eResult, tkm_tErr_SurfaceNotLoaded );
  
  /* load the vertex set and check if it was loaded. */
  DebugNote( ("Loading vertex set") );
  eSurface = Surf_LoadVertexSet( gSurface[iType], isName, iSet );
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
         eResult, tkm_tErr_CouldntLoadSurfaceVertexSet );
  
  DebugNote( ("Checking if vertex set is loaded") );
  eSurface = Surf_IsVertexSetLoaded( gSurface[iType], iSet, &bLoaded );
  DebugAssertThrowX( (bLoaded),eResult, tkm_tErr_CouldntLoadSurfaceVertexSet );
  
  /* get command and flag to set */
  switch( iSet ) {
  case Surf_tVertexSet_Pial:
    flag    = DspA_tDisplayFlag_PialSurface;
    command = tkm_tTclCommand_ShowPialSurfaceViewingOptions;
    break;
  case Surf_tVertexSet_Original:
    flag    = DspA_tDisplayFlag_OriginalSurface;
    command = tkm_tTclCommand_ShowOriginalSurfaceViewingOptions;
    break;
  default:
    DebugAssertThrowX( (1), eResult, tkm_tErr_CouldntLoadSurfaceVertexSet );
    break;
  }
  
  /* turn flag on or off and enable or disable viewing optiosn */
  DebugNote( ("Turning on vertex set in window") );
  MWin_SetDisplayFlag( gMeditWindow, -1, flag, bLoaded );
  DebugNote( ("Turning on viewing options in interface") );
  tkm_SendTclCommand( command, bLoaded?"1":"0" );
  
  /* set the surface in the window to purge the cache */
  DebugNote( ("Setting surface in window") );
  MWin_SetSurface( gMeditWindow, -1, tkm_tSurfaceType_Main, 
       gSurface[tkm_tSurfaceType_Main] );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void UnloadSurface ( tkm_tSurfaceType iType ) {
  
  tkm_tErr  eResult  = tkm_tErr_NoErr;
  Surf_tErr eSurface = Surf_tErr_NoErr;
  
  DebugAssertThrowX( (NULL != gSurface[iType]),
         eResult, tkm_tErr_InvalidParameter );
  if( !gSurface[iType] )
    return;
  
  /* free the surface. */
  eSurface = Surf_Delete( &gSurface[iType] );
  if( Surf_tErr_NoErr != eSurface ) {
    DebugPrint( ( "Surf error %d in UnloadSurface: %s\n",
      eSurface, Surf_GetErrorString( eSurface ) ) );
  }
  
  /* if this is our main surface, disable our loading options */
  if( tkm_tSurfaceType_Main == iType ) {
    tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, "0" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceViewingOptions, "0" );
    tkm_SendTclCommand(tkm_tTclCommand_ShowPialSurfaceViewingOptions,"0");
    tkm_SendTclCommand(tkm_tTclCommand_ShowOriginalSurfaceViewingOptions,"0");
    MWin_SetDisplayFlag( gMeditWindow, -1,
			 DspA_tDisplayFlag_MainSurface, FALSE );
    MWin_SetDisplayFlag( gMeditWindow, -1, 
			 DspA_tDisplayFlag_PialSurface, FALSE );
    MWin_SetDisplayFlag( gMeditWindow, -1, 
			 DspA_tDisplayFlag_OriginalSurface, FALSE );
    MWin_SetDisplayFlag( gMeditWindow, -1,
			 DspA_tDisplayFlag_InterpolateSurfaceVertices, FALSE );
  }
  
  /* update the medit window. */
  MWin_SetSurface( gMeditWindow, -1, iType, gSurface[iType] );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void WriteSurfaceValues ( tkm_tSurfaceType iType,
        char*       isFileName ) {
  
  tkm_tErr  eResult  = tkm_tErr_NoErr;
  
  DebugAssertThrowX( (NULL != gSurface[iType]),
         eResult, tkm_tErr_InvalidParameter );
  if( !gSurface[iType] )
    return;
  
  /* Write the values for this surface. */
  Surf_WriteValues( gSurface[iType], isFileName );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

// ===========================================================================



/* =========================================================== TCL WRAPPERS */

int TclCrashHard ( ClientData inClientData, 
       Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  char* pBadPtr = 0x0;
  *pBadPtr = 1;
  
  return TCL_OK;
}

int TclTranslateOverlayRegistration ( ClientData inClientData, 
              Tcl_Interp* inInterp,
              int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp,
        "wrong#args: TranslateOverlayRegisgtration distance x,y,z",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    TranslateOverlayRegisgtration ( atof(argv[1]), argv[2][0] ); 
  }
  
  return TCL_OK;
}

int TclRotateOverlayRegistration ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: RotateOverlayRegistration degrees x,y,z",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    RotateOverlayRegistration ( atof(argv[1]), argv[2][0] ); 
  }
  
  return TCL_OK;
}

int TclScaleOverlayRegistration ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp,
        "wrong#args: ScaleOverlayRegisgtration distance x,y,z",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    ScaleOverlayRegisgtration ( atof(argv[1]), argv[2][0] ); 
  }
  
  return TCL_OK;
}

int TclLoadHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadHeadPts points_file transform_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
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

int TclLoadGCA ( ClientData inClientData, Tcl_Interp* inInterp,
     int argc, char* argv[] ) {
  
  if ( argc < 3 || argc > 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadGCA volume_dir transform_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadGCA ( argv[1], argv[2] ); 
    if (argc == 4)
      {
  printf("mean filtering conditional densities %d times...\n",
         atoi(argv[3])) ;
  GCAmeanFilterConditionalDensities(gGCAVolume, atoi(argv[3])) ;
      }
  }
  
  return TCL_OK;
}
int TclLoadGCARenormalization ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadGCARenormalization renorm_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadGCARenormalization(argv[1]) ;
  }
  
  return TCL_OK;
}
int TclSaveGCA ( ClientData inClientData, Tcl_Interp* inInterp,
     int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveGCA gca_file_name",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SaveGCA ( argv[1] ); 
  }
  
  return TCL_OK;
}

int TclReadVoxelLabels ( ClientData inClientData, 
       Tcl_Interp* inInterp,
       int argc, char* argv[] )
{
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadVoxelLabels vli1 vli2",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadVLIs ( argv[1], argv[2] ); 
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

int TclThresholdVolume ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: ThresholdVolume threshold_value {0=below,1=above} new_value",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    ThresholdVolume( atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) );
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
    if( atoi(argv[1]) ) {
      FlipVolume( mri_tOrientation_Sagittal );
    }
    if( atoi(argv[2]) ) {
      FlipVolume( mri_tOrientation_Horizontal );
    }
    if( atoi(argv[3]) ) {
      FlipVolume( mri_tOrientation_Coronal );
    }
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
    switch( argv[1][0] ) {
    case 'x':
      RotateVolume( mri_tOrientation_Sagittal, atof( argv[2] ) );
      break;
    case 'y':
      RotateVolume( mri_tOrientation_Coronal, atof( argv[2] ) );
      break;
    case 'z':
      RotateVolume( mri_tOrientation_Horizontal, atof( argv[2] ) );
      break;
    default:
      Tcl_SetResult ( inInterp, 
          "No axis specified", TCL_VOLATILE );
      return TCL_ERROR;
    }
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
    LoadVolume( tkm_tVolumeType_Main, argv[1] );
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
    LoadVolume( tkm_tVolumeType_Aux, argv[1] );
    
    /* show the aux volume */
    MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_AuxVolume, 
       TRUE );
  }
  
  return TCL_OK;
}

int TclUnloadVolume ( ClientData inClientData, Tcl_Interp* inInterp,
		      int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UnloadVolume 1=aux",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    UnloadVolume( atoi(argv[1]) );
  }
  
  return TCL_OK;
}

int TclSaveVolume ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveVolume 0=main,1=aux",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SaveVolume( atoi( argv[1] ), NULL );
  }
  
  return TCL_OK;
}

int TclSaveVolumeAs ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveVolumeAs 0=main,1=aux "
		    "volume_path:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SaveVolume( atoi( argv[1] ), argv[2] );
  }
  
  return TCL_OK;
}

int TclLoadVolumeDisplayTransform ( ClientData inClientData, 
            Tcl_Interp* inInterp,
            int argc, char* argv[] ) {
  
  tkm_tVolumeType volume = tkm_tVolumeType_Main;
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadVolumeDisplayTransform "
        "volume transform:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* get a volume index */
    volume = atoi( argv[1] );
    
    /* make sure it's main or aux. if we have that volume, load the
       transform. */
    if( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if( NULL != gAnatomicalVolume[ volume ] ) {
	LoadDisplayTransform( volume, argv[2] );
      }
    }
  }
  
  return TCL_OK;
}

int TclUnloadVolumeDisplayTransform ( ClientData inClientData, 
              Tcl_Interp* inInterp,
              int argc, char* argv[] ) {
  
  tkm_tVolumeType volume = tkm_tVolumeType_Main;
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UnloadVolumeDisplayTransform "
        "volume", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* get a volume index */
    volume = atoi( argv[1] );
    
    /* make sure it's main or aux. if we have that volume, unload the
       transform. */
    if( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if( NULL != gAnatomicalVolume[ volume ] ) {
	UnloadDisplayTransform( volume );
      }
    }
  }
  
  return TCL_OK;
}
int TclUnloadGCA ( ClientData inClientData, 
       Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UnloadGCA "
        "volume", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    if (gGCAVolume)
      GCAfree(&gGCAVolume) ;
    if (gGCATransform)
      TransformFree(&gGCATransform) ;
  }
  
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

int TclClearUndoVolume ( ClientData inClientData, 
       Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ClearUndoVolume",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    ClearUndoVolume();
  }
  
  return TCL_OK;
}

int TclSetVolumeColorScale ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  tkm_tVolumeType volume = tkm_tVolumeType_Main;
  
  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeColorScale volume threshold:float squash:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* get a volume index */
    volume = atoi( argv[1] );
    
    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if( NULL != gAnatomicalVolume[ volume ] ) {
	SetVolumeBrightnessAndContrast( volume,
					atof( argv[2] ), atof( argv[3] ));
      }
    }
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

int TclSelectVoxelsByFuncValue ( ClientData inClientData, Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  FunV_tFindStatsComp compareType;
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SelectVoxelsByFuncValue "
        "compare_tpye",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    compareType = atoi( argv[1] );
    if( compareType > FunV_tFindStatsComp_Invalid &&
  compareType < FunV_knNumFindStatsComp )
      SelectVoxelsByFuncValue( compareType );
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
    GotoSavedCursor();
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
  
  xVoxel cursor;
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: NewControlPoint",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    MWin_GetCursor ( gMeditWindow, &cursor );
    NewControlPoint ( &cursor, TRUE );
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
    WriteControlPointFile ();
  }
  
  return TCL_OK;
}

int TclLoadMainSurface ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;
  
  switch( argc ) {
  case 2:
    type = tkm_tSurfaceType_Main;
    sFileName = argv[1];
    break;
  case 3:
    type = atoi( argv[1] );
    sFileName = argv[2];
    break;
  default:
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadMainSurface 0=main,1=aux"
        "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadSurface ( type, sFileName );
  }
  
  return TCL_OK;
}

int TclLoadPialSurface ( ClientData inClientData, Tcl_Interp* inInterp,
            int argc, char* argv[] ) {
  
  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;
  
  switch( argc ) {
  case 2:
    type = tkm_tSurfaceType_Main;
    sFileName = argv[1];
    break;
  case 3:
    type = atoi( argv[1] );
    sFileName = argv[2];
    break;
  default:
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadPialSurface 0=main,1=aux"
        "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( type, Surf_tVertexSet_Pial, sFileName );
  }
  
  return TCL_OK;
}

int TclLoadOriginalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;
  
  switch( argc ) {
  case 2:
    type = tkm_tSurfaceType_Main;
    sFileName = argv[1];
    break;
  case 3:
    type = atoi( argv[1] );
    sFileName = argv[2];
    break;
  default:
    Tcl_SetResult ( inInterp, 
        "wrong # args: LoadOriginalSurface 0=main,1=aux"
        "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( type, Surf_tVertexSet_Original, sFileName );
  }  
  
  return TCL_OK;
}

int TclUnloadSurface ( ClientData inClientData, Tcl_Interp* inInterp,
		       int argc, char* argv[] ) {
  
  tkm_tSurfaceType type = tkm_tSurfaceType_Main;
  
  switch( argc ) {
  case 1:
    type = tkm_tSurfaceType_Main;
    break;
  case 2:
    type = atoi( argv[1] );
    break;
  default:
    Tcl_SetResult ( inInterp, "wrong # args: UnloadSurface 0=main,1=aux",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    UnloadSurface ( type );
  }  
  
  return TCL_OK;
}

int TclUnloadAllSurfaces ( ClientData inClientData, Tcl_Interp* inInterp,
			   int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: UnloadAllSurfaces",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    UnloadSurface ( tkm_tSurfaceType_Main );
    UnloadSurface ( tkm_tSurfaceType_Aux );
  }  
  
  return TCL_OK;
}


int TclWriteSurfaceValues ( ClientData inClientData, Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  tkm_tSurfaceType type       = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;
  
  switch( argc ) {
  case 2:
    type = tkm_tSurfaceType_Main;
    sFileName = argv[1];
    break;
  case 3:
    type = atoi( argv[1] );
    sFileName = argv[2];
    break;
  default:
    Tcl_SetResult ( inInterp, 
        "wrong # args: WriteSurfaceValues 0=main,1=aux "
        "file_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    WriteSurfaceValues ( type, sFileName );
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

int TclGotoPialVertex ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: GotoPialVertex vertex_num:int",
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

int TclShowNearestPialVertex ( ClientData inClientData, 
            Tcl_Interp* inInterp,
            int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestPialVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    FindNearestSurfaceVertex ( Surf_tVertexSet_Pial );
  }  
  
  return TCL_OK;
}


int TclShowNearestInterpolatedMainVertex ( ClientData inClientData,
					   Tcl_Interp* inInterp,
					   int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestInterpolatedMainVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet_Main );
  }  
  
  return TCL_OK;
}

int TclShowNearestInterpolatedOriginalVertex ( ClientData inClientData, 
					       Tcl_Interp* inInterp,
					       int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestInterpolatedOriginalVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet_Original );
  }  
  
  return TCL_OK;
}

int TclShowNearestInterpolatedPialVertex ( ClientData inClientData, 
						Tcl_Interp* inInterp,
						int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: ShowNearestInterpolatedPialVertex",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet_Pial );
  }  
  
  return TCL_OK;
}

int TclAverageSurfaceVertexPositions ( ClientData inClientData, 
				       Tcl_Interp* inInterp,
				       int argc, char* argv[] ) {
  
  if ( argc != 2) {
    Tcl_SetResult ( inInterp, "wrong # args: AverageSurfaceVertexPositions "
		    "num_averages", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    AverageSurfaceVertexPositions ( atoi(argv[1]) );
  }  
  
  return TCL_OK;
}


int TclNewSegmentationVolume ( ClientData inClientData, 
			       Tcl_Interp* inInterp,
			       int argc, char* argv[] ) {
  
  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: NewSegmentationVolume "
		    "volume:0=main,1=aux from_anatomical:0=main,1=aux "
		    "color_file:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    NewSegmentationVolume ( atoi(argv[1]), atoi(argv[2]), argv[3] );
  }  
  
  return TCL_OK;
}

int TclLoadSegmentationVolume ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadSegmentationVolume "
		    "volume:0=main,1=aux directory_and_prefix:string "
		    "color_file:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    LoadSegmentationVolume ( atoi(argv[1]), argv[2], argv[3] );
  }  
  
  return TCL_OK;
}

int TclRecomputeSegmentation ( ClientData inClientData, 
			       Tcl_Interp* inInterp,
			       int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RecomputeSegmentation "
		    "volume:0=main,1=aux", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    RecomputeSegmentation( atoi(argv[1]) );
  }  

  return TCL_OK;
}

int TclImportSurfaceAnnotationToSegmentation ( ClientData inClientData, 
					       Tcl_Interp* inInterp,
					       int argc, char* argv[] ) {
  
  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: "
		    "ImportSurfaceAnnotationToSegmentation "
		    "volume:[main=0,aux=1] annotation:string "
		    "color_file:string", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    ImportSurfaceAnnotationToSegmentation ( atoi(argv[1]), argv[2], argv[3] );
  }  
  
  return TCL_OK;
}

int TclSetGCADisplayStatus ( ClientData inClientData, 
           Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetGCADisplayStatus new_status",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    gDisplayIntermediateResults = atoi(argv[1]) ;
  }  
  
  return TCL_OK;
}

int TclRestorePreviousSegmentation ( ClientData inClientData, 
				     Tcl_Interp* inInterp,
				     int argc, char* argv[] ) {
  
  tkm_tErr     eResult   = tkm_tErr_NoErr;
  Volm_tErr    eVolume   = Volm_tErr_NoErr;
  tkm_tSegType volume    = tkm_tSegType_Main;

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: RestorePreviousSegmentation "
		    "volume:0=main,1=aux", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  /* Check the volume parameter since we're about to use it. */
  volume = (tkm_tSegType) atoi(argv[1]);
  DebugAssertThrowX( (volume >= 0 && volume < tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );

  if( gbAcceptingTclCommands && 
      NULL != gPreviousSegmentationVolume[volume] ) {
    
    /* Delete seg volume. */
    DebugNote( ("Deleting segmentation") );
    eVolume = Volm_Delete( &gSegmentationVolume[volume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), eVolume, 
		       tkm_tErr_ErrorAccessingSegmentationVolume );
    
    /* Copy backup volume into seg volume. */
    DebugNote( ("Copying backup segmentation") );
    eVolume = Volm_DeepClone( gPreviousSegmentationVolume[volume],
			      &gSegmentationVolume[volume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), eVolume, 
		       tkm_tErr_ErrorAccessingSegmentationVolume );

    UpdateAndRedraw() ;
  }  
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  DebugCatchError( eVolume, Volm_tErr_NoErr, Volm_GetErrorString );
  tkm_DisplayError( "Recomputing Segemention", Volm_GetErrorString(eVolume),
		    "Tkmedit couldn't restore previous segmentation") ;
  EndDebugCatch;

  return TCL_OK;
}

int TclSaveSegmentationVolume ( ClientData inClientData, 
				Tcl_Interp* inInterp,
				int argc, char* argv[] ) {
  
  if ( argc != 2 && argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveSegmentationVolume "
		    "volume:0=main,1=aux [directory:string]", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {

    /* if they gave us a filename, use it, else pass null. */
    if( argc == 2 ) {
      SaveSegmentationVolume( atoi(argv[1]), NULL );
    } else {
      SaveSegmentationVolume( atoi(argv[1]), argv[2] );
    }
  }  
  
  return TCL_OK;
}

int TclExportChangedSegmentationVolume ( ClientData inClientData, 
					 Tcl_Interp* inInterp,
					 int argc, char* argv[] ) {
  
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ExportChangedSegmentationVolume "
		    "volume:0=main,1=aux [directory:string]", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    ExportChangedSegmentationVolume( atoi(argv[1]), argv[2] );
  }  
  
  return TCL_OK;
}

int TclSetSegmentationAlpha ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: SetSegmentationAlpha alpha:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SetSegmentationAlpha ( atof(argv[1]) );
  }  
  
  return TCL_OK;
}

int TclLoadFunctionalOverlay ( ClientData inClientData, 
             Tcl_Interp* inInterp,
             int argc, char* argv[] ) {
  
  char sPathAndStem[tkm_knPathLen] = "";
  char sRegistration[tkm_knPathLen] = "";
  
  if ( argc != 2 && argc != 3 && argc != 4 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: LoadFunctionalOverlay directory "
        "stem [registration]",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* if we have a path and a stem, first concatenate the first two
       arguments to make a stem and path.  then, if we have a third
       argument, pass it as the registration, else pass null. */
    if( argc >= 3 && argv[2][0] != '\0' ) {
      xUtil_snprintf( sPathAndStem, sizeof(sPathAndStem),
		      "%s/%s", argv[1], argv[2] );
    } else {
      /* no path and stem, just use the file name. */
      xUtil_strncpy( sPathAndStem, argv[1], sizeof(sPathAndStem) );
    }

    if( argc == 4 &&
	argv[3][0] != '\0' ) {
      xUtil_strncpy( sRegistration, argv[3], sizeof(sRegistration) );
      LoadFunctionalOverlay( sPathAndStem, NULL, sRegistration );
    } else {
      LoadFunctionalOverlay( sPathAndStem, NULL, NULL );
    }
  }  
  
  return TCL_OK;
}

int TclLoadFunctionalTimeCourse ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  char sPathAndStem[tkm_knPathLen]  = "";
  char sRegistration[tkm_knPathLen] = "";
  
  if ( argc != 3 && argc != 4 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: LoadFunctionalTimeCourse directory "
        "stem [registration]",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* first concatenate the first two arguments to make a stem and path.
       then, if we have a third argument, pass it as the registration,
       else pass null. */
    xUtil_snprintf( sPathAndStem,  sizeof(sPathAndStem),
        "%s/%s", argv[1], argv[2] );
    if( argc == 4 &&
  argv[3][0] != '\0' ) {
      xUtil_strncpy( sRegistration, argv[3], sizeof(sRegistration) );
      LoadFunctionalTimeCourse( sPathAndStem, NULL, sRegistration );
    } else {
      LoadFunctionalTimeCourse( sPathAndStem, NULL, NULL );
    }
  }  
  
  return TCL_OK;
}

int TclLoadDTIVolumes ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  if ( argc != 6 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadDTIVolumes "
		    "ev_volume_name:string fa_volume_name:string "
                    "red_axis:[0=x,y=1,z=2] "
		    "green_axis:[0=x,y=1,z=2] blue_axis:[0=x,y=1,z=2]",
		    TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    eResult = LoadDTIVolume( argv[1], argv[2],
			     atoi(argv[3]), atoi(argv[4]), atoi(argv[5]) );
  }
  
  return TCL_OK;
}

int TclSetDTIAlpha ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: SetDTIAlpha alpha:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SetDTIAlpha ( atof(argv[1]) );
  }  
  
  return TCL_OK;
}

int TclSmoothFunctionalOverlay ( ClientData inClientData, 
         Tcl_Interp* inInterp,
         int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: SmoothFunctionalOverlay sigma",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    SmoothOverlayData( atof(argv[1]) );
  }  
  
  return TCL_OK;
}

int TclSetTimerStatus ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {
  
  tBoolean bStatus = FALSE;
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, 
        "wrong # args: SetTimerStatus 0|1", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    bStatus = (tBoolean) atoi( argv[1] );
    if( bStatus && !gbTimerOn ) {
      StartTimer();
    } else if( !bStatus && gbTimerOn ) {
      StopTimer();
    }
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

int TclRedrawScreen ( ClientData inClientData, 
          Tcl_Interp* inInterp,
          int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: RedrawScreen",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    if( gMeditWindow != NULL ) {
      MWin_ForceRedraw( gMeditWindow );
    }
  }  
  
  return TCL_OK;
}

int TclDebugPrint ( ClientData inClientData, 
        Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp,
        "wrong # args: DebugPrint string",
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  DebugPrint( ( "Tcl: %s\n", argv[1] ) );
  
  return TCL_OK;
}

int TclStartTimer ( ClientData inClientData, Tcl_Interp* inInterp,
        int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: StartTimer", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  StartTimer();
  
  return TCL_OK;
}

int TclStopTimer ( ClientData inClientData, Tcl_Interp* inInterp,
       int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: StopTimer", TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  StopTimer();
  
  return TCL_OK;
}

int TclExecuteQueuedScripts ( ClientData inClientData, Tcl_Interp* inInterp,
            int argc, char* argv[] ) {
  
  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ExecuteQueuedScripts", 
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  ExecuteQueuedScripts();
  
  return TCL_OK;
}

int TclGetSubjectName ( ClientData inClientData, Tcl_Interp* inInterp,
      int argc, char* argv[] ) {
  
  Volm_tErr    eVolume          = tkm_tErr_NoErr;
  char      sSubjectName[tkm_knPathLen] = "";
  tkm_tVolumeType volume          = tkm_tVolumeType_Main;
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: GetSubjectName 0=main|1=aux", 
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* get the volume type and check it */
    volume = (tkm_tVolumeType) atoi(argv[1]);
    if( volume != tkm_tVolumeType_Main && volume != tkm_tVolumeType_Aux ) {
      Tcl_SetResult ( inInterp, "wrong # args: GetSubjectName 0=main|1=aux", 
          TCL_VOLATILE );
      return TCL_ERROR;
    }
    
    /* get the subject name */
    eVolume = Volm_CopySubjectName( gAnatomicalVolume[volume], sSubjectName,
            sizeof( sSubjectName ));
    if( Volm_tErr_NoErr != eVolume ) {
      Tcl_SetResult ( inInterp, "Couldn't get volume name. Is volume loaded?",
          TCL_VOLATILE );
      return TCL_ERROR;
    }
    
    /* return it */
    Tcl_SetResult( inInterp, sSubjectName, TCL_VOLATILE );
  }
  
  return TCL_OK;
}

int TclGetSubjectDir ( ClientData inClientData, Tcl_Interp* inInterp,
           int argc, char* argv[] ) {
  
  Volm_tErr    eVolume          = tkm_tErr_NoErr;
  char      sSubjectDir[tkm_knPathLen]  = "";
  tkm_tVolumeType volume          = tkm_tVolumeType_Main;
  
  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: GetSubjectDir 0=main|1=aux", 
        TCL_VOLATILE );
    return TCL_ERROR;
  }
  
  if( gbAcceptingTclCommands ) {
    
    /* get the volume type and check it */
    volume = (tkm_tVolumeType) atoi(argv[1]);
    if( volume != tkm_tVolumeType_Main && volume != tkm_tVolumeType_Aux ) {
      Tcl_SetResult ( inInterp, "wrong # args: GetSubjectDir 0=main|1=aux", 
          TCL_VOLATILE );
      return TCL_ERROR;
    }
    
    /* get the source dir */
    eVolume = Volm_CopySourceDir( gAnatomicalVolume[volume], sSubjectDir,
          sizeof( sSubjectDir ));
    if( Volm_tErr_NoErr != eVolume ) {
      Tcl_SetResult ( inInterp, "Couldn't get volume directory. Is volume loaded?", TCL_VOLATILE );
      return TCL_ERROR;
    }
    
    /* return it */
    Tcl_SetResult( inInterp, sSubjectDir, TCL_VOLATILE );
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

static char sTclEnvVar[STRLEN] = "";
static char sTkEnvVar[STRLEN] = "";


int main ( int argc, char** argv ) {
  
  tkm_tErr   eResult                           = tkm_tErr_NoErr;
  x3Lst_tErr e3DList                           = x3Lst_tErr_NoErr;
  FunV_tErr  eFunctional                       = FunV_tErr_NoError;
  MWin_tErr  eWindow                           = MWin_tErr_NoErr;
  tkm_tVolumeType volume                       = tkm_tVolumeType_Main;
  int        eTcl                              = TCL_OK;
  tBoolean   bFoundInterface                   = FALSE;
  char       sInterfaceFileName[tkm_knPathLen] = "";
  char*      pEnvVar                           = NULL;
  FILE*      pFile                             = NULL ;
  int        nArg                              = 0;
  time_t     theTime;
  char       sSubjectName[tkm_knNameLen]       = "";
  
#ifdef SET_TCL_ENV_VAR
  tBoolean  bChangedEnvVar    = FALSE;
  char*      sTclLib        = NULL;
  char      sSavedTclLib[STRLEN] = "";
  int      nTclLength        = 0;
  char      sNewTclLib[STRLEN]   = "";
  char*      sTkLib        = NULL;
  char      sSavedTkLib[STRLEN]   = "";
  int      nTkLength        = 0;
  char      sNewTkLib[STRLEN]   = "";
#endif
  
  /* init our debugging macro code, if any. */
  InitDebugging( "tkmedit" );
  EnableDebuggingOutput;
  
  /* install our segfault handler */
  DebugRegisterSegfaultHandler( HandleSegfault );
  
  //  xDbg_PrintStatus();
  
  DebugEnterFunction( ("main()") );
  
  gm_screen2ras = MatrixAlloc(4,4,MATRIX_REAL) ;
  *MATRIX_RELT(gm_screen2ras, 1, 1) = -1 ;
  *MATRIX_RELT(gm_screen2ras, 2, 3) = 1 ;
  *MATRIX_RELT(gm_screen2ras, 3, 2) = -1 ;
  
  *MATRIX_RELT(gm_screen2ras, 1, 4) = 128 ;
  *MATRIX_RELT(gm_screen2ras, 2, 4) = -128 ;
  *MATRIX_RELT(gm_screen2ras, 3, 4) = 128 ;
  
  *MATRIX_RELT(gm_screen2ras, 4, 4) = 1 ;
  
  
  gm_ras2screen = MatrixInverse(gm_screen2ras, NULL) ;

  /* write the time started, progam name, and arguments to the debug output */
  time( &theTime );
  DebugPrint( ( "tkmedit started: %s\n\t", ctime( &theTime ) ) );
  for( nArg = 0; nArg < argc; nArg++ ) {
    DebugPrint( ( "%s ", argv[nArg] ) );
  }
  DebugPrint( ( "\n\n" ) );
  
  /* init glut */
  DebugNote( ("Initializing glut") );
  glutInit( &argc, argv );
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
  
  /* init main window */
  DebugNote( ("Creating main window") );
  eWindow = MWin_New( &gMeditWindow, "", 512, 512 );
  DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );
  
  /* init our control pt list */
  DebugNote( ("Initializing control point list") );
  e3DList = x3Lst_New( &gControlPointList, NUM_UNDOS );
  DebugAssertThrow( (x3Lst_tErr_NoErr == e3DList) );
  x3Lst_SetComparator( gControlPointList, CompareVoxels );
  
  /* init other things */
  DebugNote( ("Initalizing undo list") );
  eResult = InitUndoList();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );
  DebugNote( ("Initalizing selection module") );
  eResult = InitSelectionModule();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );
  
  /* create functional volume */
  DebugNote( ("Creating functional volume") );
  eFunctional = FunV_New( &gFunctionalVolume,
        UpdateAndRedraw, tkm_SendTclCommand, 
        SendTCLCommand );
  DebugAssertThrow( (FunV_tErr_NoError == eFunctional) );
  
  /* set windows data sources */
  DebugNote( ("Setting control points space in window") );
  eWindow = MWin_SetControlPointsSpace( gMeditWindow, -1, gControlPointList );
  DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );
  
  DebugNote( ("Setting selection list in window.") );
  eWindow = MWin_SetSelectionSpace( gMeditWindow, -1, gSelectedVoxels );
  DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );
  
  /* start by disabling a bunch of stuff. if it gets loaded, the options
     will be enabled later */
  DebugNote( ("Initializing tcl display state") );
  tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceLoadingOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSurfaceViewingOptions,"0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowPialSurfaceViewingOptions,"0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowOriginalSurfaceViewingOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowDTIOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowVLIOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowGCAOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowSegmentationOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowAuxSegmentationOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowVolumeDirtyOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeDirtyOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowMainTransformLoadedOptions, "0" );
  tkm_SendTclCommand( tkm_tTclCommand_ShowAuxTransformLoadedOptions, "0" );
  
  /* init variables */
  for( volume = 0; volume < tkm_knNumVolumeTypes; volume++ ) {
    gAnatomicalVolume[volume] = NULL;
    SetVolumeDirty( volume, FALSE );
  }

  for( volume = 0; volume < tkm_knNumSegTypes; volume++ ) {
    gPreviousSegmentationVolume[volume] = NULL;
    gSegmentationVolume[volume] = NULL;
    gSegmentationChangedVolume[volume] = NULL;
  }
  
  /* find the user's home directory, or where the program is running from */
  DebugNote( ("Finding home directory") );
  eResult = FindUserHomeDir();
  DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
  
  /* start program, now as function; gl window not opened yet */
  DebugNote( ("Parsing command line arguments") );
  ParseCmdLineArgs( argc, argv );
  
#ifdef USE_LICENSE
  checkLicense(envptr);
#endif
  
  
  /* if there is no interface name, i.e. if they didn't pass one in on
     the command line, look in the local directory for the script, then
     in MRI_DIR/lib/tcl. if we can't fine, we have to exit. */
  pFile = NULL;
  bFoundInterface = FALSE;

  /* did they request one? */
  if( !MATCH( gInterfaceScriptName, "" ) ) {

    /* try to open it. */
    xUtil_strncpy( sInterfaceFileName, gInterfaceScriptName, 
		   sizeof(sInterfaceFileName) );
    DebugNote( ( "Trying to open %s\n", sInterfaceFileName ) );
    pFile = fopen ( sInterfaceFileName,"r" );
    if( NULL != pFile ) {
      bFoundInterface = TRUE;
	fclose( pFile );
    } else {
      tkm_DisplayError( "Opening specified interface file",
			"No valid file found",
			"Tkmedit couldn't find or open the interface file "
			"you specified. It will look elsewhere for the "
			"standard one." );
    }
  }
    
  if( !bFoundInterface ) {

    /* next, look in TKMEDIT_SCRIPTS_DIR */
    DebugNote( ("Getting TKMEDIT_SCRIPTS_DIR env var") );
    pEnvVar = getenv("TKMEDIT_SCRIPTS_DIR");
    if( NULL != pEnvVar) {
      xUtil_snprintf( sInterfaceFileName, sizeof(sInterfaceFileName),
		      "%s/tkmedit.tcl", pEnvVar); 
      DebugNote( ( "Trying to open %s", sInterfaceFileName ) );
      pFile = fopen( sInterfaceFileName,"r" );
      if( NULL != pFile ) {
	bFoundInterface = TRUE;
	fclose( pFile );
      } else {
	tkm_DisplayError( "Opening interface file in TKMEDIT_SCRIPTS_DIR",
			  "No valid file found",
			  "Tkmedit couldn't find or open the interface file "
			  "in TKMEDIT_SCRIPTS_DIR. It will look elsewhere "
			  "for the standard one." );
      }
    }
  }    
  
  
  if( !bFoundInterface ) {

    /* look in the local directory */
    DebugNote( ("Trying local tkmedit.tcl") );
    xUtil_strncpy( sInterfaceFileName, "tkmedit.tcl", 
		   sizeof(sInterfaceFileName) ); 
    pFile = fopen( sInterfaceFileName, "r" );
    if( NULL != pFile ) {
      bFoundInterface = TRUE;
      fclose( pFile );
    } 
  }

  if ( !bFoundInterface ) {

    /* finally try in MRI_DIR/lib/tcl. make sure we have MRI_DIR
       defined */
    DebugNote( ("Getting MRI_DIR env var") );
    pEnvVar = getenv("MRI_DIR");
    if( NULL == pEnvVar) {
      tkm_DisplayError( "Trying to find interface file",
			"No valid file found",
			"Tkmedit couldn't find a valid interface file. "
			"Normally this is in the directory specified in "
			"the MRI_DIR varible, but this was not set in "
			"your environment. Tkmedit needs this file to "
			"run." );
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }
    
    xUtil_snprintf( sInterfaceFileName, sizeof(sInterfaceFileName),
		    "%s/lib/tcl/%s", pEnvVar, "tkmedit.tcl"); 
    DebugNote( ( "Trying to open %s", sInterfaceFileName ) );
    pFile = fopen( sInterfaceFileName,"r" );
    if( NULL != pFile ) {
      bFoundInterface = TRUE;
      fclose( pFile );
    } 
  }

  
  /* if file still not found bail out. */
  if ( !bFoundInterface ) {
    tkm_DisplayError( "Trying to find interface file",
          "No valid file found",
          "Tkmedit couldn't find a valid interface file. "
          "Normally this is in the directory specified in "
          "the MRI_DIR varible, but it can also be in the same "
          "directory as tkmedit. Tkmedit needs this file to "
          "run. Try reinstalling your distribution." );
    PrintCachedTclErrorDlogsToShell();
    exit( 1 );
  }
  
  DebugPrint( ( "Using interface file %s\n", sInterfaceFileName) );
  
  /* process ctrl pts file */
  DebugNote( ("Processing control point file.") );
  ProcessControlPointFile();
  
  /* set window's data */
  DebugNote( ("Setting window title.") );
  Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main], 
      sSubjectName, sizeof(sSubjectName) );
  MWin_SetWindowTitle( gMeditWindow, sSubjectName );
  
  DebugNote( ("Initalizing undo volume") );
  eResult = InitUndoVolume();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );
  
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
	DebugPrint( ( "Couldn't set TCL_LIBRARY to %s\n", sNewTclLib
		      ) );
	exit( 1 );
      }
      sprintf( sTkEnvVar, "%s=%s", "TK_LIBRARY", sNewTkLib );
      if( putenv( sTkEnvVar ) ) {
  OutputPrint "ERROR: Couldn't set TK_LIBRARY env var.\n" 
    EndOutputPrint;
  DebugPrint( ( "Couldn't set TK_LIBRARY to %s\n", sNewTkLib
          ) );
  exit( 1 );
      }
    }
    
  } else {
    
    OutputPrint "ERROR: TCL_LIBRARY or TK_LIBRARY environement variable is not set.\n" EndOutputPrint;
    DebugPrint( ( "TCL_LIBRARY or TK_LIBRARY env var not set.\n" ) );exit ( 1 );
  }
  
#endif
  
  /* start tcl/tk; first make interpreter */
  DebugNote( ("Creating Tcl interpreter") );
  interp = Tcl_CreateInterp();
  
  // kt - set global interp
  DebugNote( ("Setting global interpreter") );
  SetTclInterp ( interp );
  
  /* read tcl/tk internal startup scripts */
  eTcl = Tcl_Init( interp );
  if( TCL_OK != eTcl ) {
    DebugPrint( ("Tcl_Init returned %d\n", (int)eTcl) );
    tkm_DisplayError( "Initializing Tcl",
          "Error initializing Tcl",
          "For some reason, Tcl couldn't be initialized. Possible "
          "reasons include it not being installed, installed "
          "incorrectly, or the TCL_LIBRARY environment variable "
          "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTcl );
  }
  eTcl = Tk_Init( interp );
  if( TCL_OK != eTcl ) {
    DebugPrint( ("Tk_Init returned %d\n", (int)eTcl) );
    tkm_DisplayError( "Initializing Tk",
          "Error initializing Tk",
          "For some reason, Tk couldn't be initialized. Possible "
          "reasons include it not being installed, installed "
          "incorrectly, or the TK_LIBRARY environment variable "
          "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTk );
  }
  eTcl = Tix_Init( interp );
  if( TCL_OK != eTcl ) {
    DebugPrint( ("Tix_Init returned %d\n", (int)eTcl) );
    tkm_DisplayError( "Initializing Tix",
          "Error initializing Tix",
          "For some reason, Tix couldn't be initialized. Possible "
          "reasons include it not being installed, installed "
          "incorrectly, or the TIX_LIBRARY environment variable "
          "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTix );
  }
  eTcl = Blt_Init( interp );
  if( TCL_OK != eTcl ) {
    DebugPrint( ("Blt_Init returned %d\n", (int)eTcl) );
    tkm_DisplayError( "Initializing BLT",
          "Error initializing BLT",
          "For some reason, BLT couldn't be initialized. Possible "
          "reasons include it not being installed, installed "
          "incorrectly, or the TIX_LIBRARY environment variable "
          "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitBLT );
  }
  
  Tcl_StaticPackage( interp, "BLT", Blt_Init, Blt_SafeInit );
  Tcl_StaticPackage( interp, "Tix", Tix_Init, Tix_SafeInit );
  
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
  DebugNote( ("Registereing window tcl commands.") );
  MWin_RegisterTclCommands ( gMeditWindow, interp );
  
  /* if we're running in a terminal shell, make our tcl shell interactive
     so we can enter commands. */
  DebugNote( ("Determining if this is a tty shell") );
  tty = isatty(0);
  DebugNote( ("Setting tcl_interactive var to %d", (int)tty) );
  Tcl_SetVar( interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY );
  
  /* register tcl commands */
  DebugNote( ("Registering tkmedit tcl commands") );
  
  Tcl_CreateCommand ( interp, "SetGCADisplayStatus", TclSetGCADisplayStatus,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "RecomputeSegmentation", TclRecomputeSegmentation,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "RestorePreviousSegmentation", 
          TclRestorePreviousSegmentation,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "CrashHard", TclCrashHard,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ReadVoxelLabels", TclReadVoxelLabels,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "DebugPrint", TclDebugPrint,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RotateOverlayRegistration",
          TclRotateOverlayRegistration,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "TranslateOverlayRegistration",
          TclTranslateOverlayRegistration,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ScaleOverlayRegistration",
          TclScaleOverlayRegistration,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
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
  
  Tcl_CreateCommand ( interp, "LoadGCA",
          TclLoadGCA,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadGCARenorm",
          TclLoadGCARenormalization,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveGCA",
          TclSaveGCA,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveRGB",
          TclSaveRGB,
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
  
  Tcl_CreateCommand ( interp, "UnloadVolume",
		      TclUnloadVolume,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveVolume",
          TclSaveVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveVolumeAs",
          TclSaveVolumeAs,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadVolumeDisplayTransform",
          TclLoadVolumeDisplayTransform,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnloadVolumeDisplayTransform",
          TclUnloadVolumeDisplayTransform,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnloadGCA",
          TclUnloadGCA,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SnapshotVolume",
          TclSnapshotVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RestoreVolumeFromSnapshot",
          TclRestoreVolumeFromSnapshot,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ClearUndoVolume",
          TclClearUndoVolume,
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
  
  Tcl_CreateCommand ( interp, "SelectVoxelsByFuncValue",
          TclSelectVoxelsByFuncValue,
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
  
  Tcl_CreateCommand ( interp, "WriteControlPointFile",
          TclWriteControlPointFile,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadMainSurface",
          TclLoadMainSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadPialSurface",
          TclLoadPialSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadOriginalSurface",
          TclLoadOriginalSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnloadSurface",
          TclUnloadSurface,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "UnloadAllSurfaces",
          TclUnloadAllSurfaces,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "WriteSurfaceValues",
          TclWriteSurfaceValues,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GotoMainVertex",
          TclGotoMainVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "GotoPialVertex",
          TclGotoPialVertex,
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
  
  Tcl_CreateCommand ( interp, "ShowNearestPialVertex",
          TclShowNearestPialVertex,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedMainVertex",
		      TclShowNearestInterpolatedMainVertex,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedOriginalVertex",
		      TclShowNearestInterpolatedOriginalVertex,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedPialVertex",
		      TclShowNearestInterpolatedPialVertex,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "AverageSurfaceVertexPositions",
		      TclAverageSurfaceVertexPositions,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "NewSegmentationVolume",
		      TclNewSegmentationVolume,
		      (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "LoadSegmentationVolume",
          TclLoadSegmentationVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SaveSegmentationVolume",
          TclSaveSegmentationVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ExportChangedSegmentationVolume",
          TclExportChangedSegmentationVolume,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "ImportSurfaceAnnotationToSegmentation",
          TclImportSurfaceAnnotationToSegmentation,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SetSegmentationAlpha",
          TclSetSegmentationAlpha,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadFunctionalOverlay",
          TclLoadFunctionalOverlay,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadFunctionalTimeCourse",
          TclLoadFunctionalTimeCourse,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "LoadDTIVolumes",
          TclLoadDTIVolumes,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SetDTIAlpha",
          TclSetDTIAlpha,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SmoothFunctionalOverlay",
          TclSmoothFunctionalOverlay,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "SetTimerStatus",
          TclSetTimerStatus,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "QuitMedit",
          TclQuitMedit,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand ( interp, "RedrawScreen",
          TclRedrawScreen,
          (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand( interp, "StartTimer", TclStartTimer,
         (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand( interp, "StopTimer", TclStopTimer,
         (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand( interp, "ExecuteQueuedScripts", TclExecuteQueuedScripts,
         (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  Tcl_CreateCommand( interp, "GetSubjectName", TclGetSubjectName,
         (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand( interp, "GetSubjectDir", TclGetSubjectDir,
         (ClientData) NULL, (Tcl_CmdDeleteProc*) NULL );
  
  /* parse the interface file */
  DebugNote( ("Parsing %s", sInterfaceFileName) );
  eTcl = Tcl_EvalFile( interp, sInterfaceFileName );
  if( *interp->result != 0 ) {
    DebugPrint( ( "Error reading %s\n\t%s\n", 
      sInterfaceFileName, interp->result ) );
    tkm_DisplayError( "Parsing interface file",
          "Tcl error while parsing file",
          "There was an error in the interface script file and "
          "it could not be parsed. Unless you have modified the "
          "file, this is not your fault. Please report this "
          "error and try reinstalling the file." );
    exit( 1 );
  }
  
  /* always start up command line shell too */
  DebugNote( ("Creating TK file handler") );
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);
  
  /* finish initializing functional volume */
  DebugNote( ("Initializing graph window") );
  eFunctional = FunV_InitGraphWindow( gFunctionalVolume, interp );
  if( FunV_tErr_NoError != eFunctional ) {
    exit( 1 );
  }
  
  DebugNote( ("Registering functional Tcl commands") );
  eFunctional = FunV_RegisterTclCommands( gFunctionalVolume, interp );
  if( FunV_tErr_NoError != eFunctional ) {
    exit( 1 );
  }
  
  /* set data in window */
  DebugNote( ("Setting functional volume in main window") );
  MWin_SetOverlayVolume( gMeditWindow, -1, gFunctionalVolume );
  
#if 0
  /* show the tal coords and hide the ras coords */
  if (NULL != gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues->linear_transform)
    {
      DebugNote( ("Showing Tal coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowTalCoords, "1" );
      DebugNote( ("Showing coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowRASCoords, "0" );
    }
#endif
  
  /* show the volume coords */
  DebugNote( ("Showing volume coords") );
  tkm_SendTclCommand( tkm_tTclCommand_ShowVolumeCoords, "1" );
  
  /* now let the window accept tcl commands. */
  DebugNote( ("Telling window to accept tcl commands") ); 
  gbAcceptingTclCommands = TRUE;
  MWin_AcceptTclCommands( gMeditWindow, TRUE );
  
  DebugNote( ("Telling interface to finish building") ); 
  tkm_SendTclCommand( tkm_tTclCommand_FinishBuildingInterface, "" );
  
  /* we probably get sent a few tcl commands before now, and we cached
     them, so send them now. */
  DebugNote( ( "Sending cached tcl commands" ) );
  SendCachedTclCommands ();
  
  /* send the brightness and contrast to the tcl window. */
  SendBrightnessAndContrastUpdate( tkm_tVolumeType_Main );
  if( NULL != gAnatomicalVolume[tkm_tVolumeType_Aux] ) {
    SendBrightnessAndContrastUpdate( tkm_tVolumeType_Aux );
  }
  
  /* set default segmentation alpha */
  DebugNote( ("Setting default segmentation alpha") );
  SetSegmentationAlpha( gfSegmentationAlpha );
  
  /* set default DTI alpha */
  DebugNote( ("Setting default DTI alpha") );
  SetDTIAlpha( gfDTIAlpha );
  
  /* if using csurf interface, call the func that hides a bunch of stuff. */
  if( gbUseCsurfInterface ) {
    DebugNote( ("Enabling csurf version of interface") ); 
    tkm_SendTclCommand( tkm_tTclCommand_CsurfInterface, "" );
  }
  
  /* if the tktimer autostart environment variable is declared,
     automatically start the timer here. */
  if( getenv( ksTkTimerAutoStartEnvVar ) ) {
    StartTimer();
  }
  
  /* never returns */
  DebugNote( ("Entering main loop") );
  glutMainLoop ();
  
  DebugCatch;
  DebugCatchError( eWindow, MWin_tErr_NoErr, MWin_GetErrorString );
  DebugCatchError( e3DList, x3Lst_tErr_NoErr, x3Lst_GetErrorString );
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return 0;
}

void tkm_Quit () {
  
  DebugEnterFunction( ("tkm_Quit()") );
  
  /* if the tktimer autostart environment variable is declared,
     automatically stop the timer here. */
  if( getenv( ksTkTimerAutoStartEnvVar ) ) {
    StopTimer();
  }
  /* delete window */
  DebugNote( ( "Deleting main window.\n" ) );
  MWin_Delete( &gMeditWindow );
  
  /* delete everything we allocated before */
  if( NULL != gGCAVolume ) {
    DebugNote( ("Deleting GCA.") );
    DeleteGCA();
  }
  if( NULL != gVLI1 || NULL != gVLI2 ) {
    DebugNote( ("Deleting VLIs.") );
    DeleteVLIs();
  }
  DebugNote( ("Deleteing functional volume.") );
  FunV_Delete( &gFunctionalVolume );
  DebugNote( ("Deleting selection module.") );
  DeleteSelectionModule ();
  DebugNote( ("Deleting control point list.") );
  x3Lst_Delete( &gControlPointList );
  DebugNote( ("Deleting undo volume.") );
  DeleteUndoVolume ();
  DebugNote( ("Deleting undo list.") );
  DeleteUndoList ();
  
  DebugPrint( ("Program terminated successfully.\n") );
  DebugNote( ( "Deleting debugging and exiting.\n" ) );
  DeleteDebugging;
  
  /* shut down tcl. this also quits the app. */
  SendTCLCommand ( "exit" );
  
  exit( 0 );
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
           "\n     (script that generates prompt)");
      fprintf(stderr, "%s\n", interp->result);
      goto defaultPrompt;
    }
  }
  fflush(stdout);
}

/* ================================================ Control point utilities */

/* reads the control.dat file, 
   transforms all pts from RAS space 
   to voxel space, and adds them as 
   control pts */
void ProcessControlPointFile ( ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen] = "";
  FILE*      file         = NULL;
  float      rasX         = 0;
  float      rasY         = 0;
  float      rasZ         = 0;
  int      nNumPointsRead       = 0;
  xVoxel    ras;
  xVoxel    idx;
  
  DebugEnterFunction( ("ProcessControlPointFile ()") );
  
  /* don't parse the file if we already have. */
  if( TRUE == gbParsedControlPointFile )
    DebugGotoCleanup;
  
  /* make the file name */
  DebugNote( ("Making file name from control.dat") );
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, 
    sFileName, sizeof(sFileName) );
  
  /* open for reading. position file ptr at beginning of file. */
  DebugNote( ("Opening control point file") );
  file = fopen( sFileName, "r" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );  
  
  /* while we have points left */
  while ( !feof(file) ) {
    
    /* read in some numbers */
    DebugNote( ("Reading three floats from control.dat") );
    nNumPointsRead = fscanf( file, "%f %f %f", &rasX, &rasY, &rasZ );
    DebugAssertThrowX( (nNumPointsRead == 3 || nNumPointsRead == EOF), 
           eResult, tkm_tErr_ErrorAccessingFile ); 
    
    if( EOF != nNumPointsRead ) {
      
      /* transform from ras to voxel */
      xVoxl_SetFloat( &ras, rasX, rasY, rasZ );
      eVolume = Volm_ConvertRASToIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
              &ras, &idx );
      DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
       eResult, tkm_tErr_ErrorAccessingVolume );
      
      /* add it to our cntrl points list */
      NewControlPoint( &idx, FALSE );
    }
  }    
  
  /* mark that we have processed the file, and shouldn't do it again. */
  gbParsedControlPointFile = TRUE;
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != file ) {
    DebugNote( ("Closing control point file") );
    fclose( file );
  }
  
  DebugExitFunction;
}

/* writes all control points to the
   control.dat file in RAS space */
void WriteControlPointFile ( ) {
  
  tkm_tErr   eResult          = tkm_tErr_NoErr;
  Volm_tErr  eVolume          = Volm_tErr_NoErr;
  x3Lst_tErr e3DList          = x3Lst_tErr_NoErr;
  xList_tErr eList          = xList_tErr_NoErr;
  char       sFileName[tkm_knPathLen] = "";
  FILE*       file          = NULL;
  int       nPlane          = 0;
  xListRef   list          = NULL;
  xVoxelRef  idx          = NULL;
  xVoxel     ras;
  
  /* make the file name */
  DebugNote( ("Making file name from control.dat") );
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints, 
    sFileName, sizeof(sFileName) );
  
  /* open for writing. position file ptr at beginning of file. */
  DebugNote( ("Opening control point file") );
  file = fopen( sFileName, "w" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );  
  
  /* get the ctrl pts in the list... */
  for( nPlane = 0; nPlane < gnAnatomicalDimensionZ; nPlane++ ) {
    
    /* get the list for this x value. */
    e3DList = x3Lst_GetItemsInXPlane( gControlPointList, nPlane, &list );
    DebugAssertThrowX( (e3DList == x3Lst_tErr_NoErr),
           eResult, tkm_tErr_ErrorAccessingList );
    
    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&idx )) 
     != xList_tErr_EndOfList ) {
      
      if( idx ) {
  
  /* transform to ras space. */
  eVolume = Volm_ConvertIdxToRAS(gAnatomicalVolume[tkm_tVolumeType_Main],
               idx, &ras );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* write to the file */
  fprintf( file, "%f %f %f\n", xVoxl_ExpandFloat(&ras) );
      }
    }
    
    DebugAssertThrowX( (eList == xList_tErr_EndOfList),
           eResult, tkm_tErr_ErrorAccessingList );
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != file ) {
    DebugNote( ("Closing control point file") );
    fclose( file );
  }
  
  DebugExitFunction;
}


float FindNearestControlPoint ( xVoxelRef   inVolumeVox, 
        mri_tOrientation inPlane,
        xVoxelRef*   outCtrlPt ) {
  
  tBoolean     bFound    = FALSE;
  float         fDistance  = 0;
  float         fClosestDistance = 0;
  float         fFoundDistance  = 0;
  xListRef     list    = NULL;
  xVoxelRef    voxel    = NULL;
  xVoxelRef    closestVoxel  = NULL;
  x3Lst_tErr   e3DList    = x3Lst_tErr_NoErr;
  xList_tErr   eList    = xList_tErr_NoErr;
  
  /* get the list to search in */
  switch ( inPlane ) {
  case mri_tOrientation_Coronal: 
    e3DList = x3Lst_GetItemsInZPlane( gControlPointList, 
              xVoxl_GetZ(inVolumeVox), &list );
    break;
  case mri_tOrientation_Horizontal: 
    e3DList = x3Lst_GetItemsInYPlane( gControlPointList, 
              xVoxl_GetY(inVolumeVox), &list );
    break;
  case mri_tOrientation_Sagittal: 
    e3DList = x3Lst_GetItemsInXPlane( gControlPointList, 
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
    fClosestDistance = gnAnatomicalDimensionX * gnAnatomicalDimensionY *
      gnAnatomicalDimensionZ ;
    bFound = FALSE;
    
    /* traverse the list */
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&voxel )) 
     != xList_tErr_EndOfList ) {
      
      if( voxel ) {
  
  /* get the distance to the clicked voxel... */
  fDistance = sqrt(
       ((xVoxl_GetX(inVolumeVox) - xVoxl_GetX(voxel)) * 
        (xVoxl_GetX(inVolumeVox) - xVoxl_GetX(voxel))) +
       ((xVoxl_GetY(inVolumeVox) - xVoxl_GetY(voxel)) * 
        (xVoxl_GetY(inVolumeVox) - xVoxl_GetY(voxel))) +
       ((xVoxl_GetZ(inVolumeVox) - xVoxl_GetZ(voxel)) * 
        (xVoxl_GetZ(inVolumeVox) - xVoxl_GetZ(voxel))) );
  
  /* if it's less than our max, mark the distance and copy the vox */
  if ( fDistance < fClosestDistance ) {
    fClosestDistance = fDistance;
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
      fFoundDistance = fClosestDistance;
    }
  }
  
  goto cleanup;
  
 error:
  
  if( eList != xList_tErr_NoErr )
    DebugPrint( ( "xList error %d in FindNearestControlPoint.\n", eList ) );
  
  if( e3DList != x3Lst_tErr_NoErr )
    DebugPrint( ( "x3Lst error %d in FindNearestControlPoint.\n", eList ) );
  
 cleanup:
  
  return fFoundDistance;
}

void NewControlPoint ( xVoxelRef iCtrlPt,
           tBoolean  ibWriteToFile ) {
  
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  ctrlPt  = NULL;
  
  /* allocate a copy of the voxel */
  xVoxl_New( &ctrlPt );
  xVoxl_Copy( ctrlPt, iCtrlPt );
  
  /* add the voxel to the ctrl pt space */
  e3DList = x3Lst_AddItem( gControlPointList, ctrlPt, ctrlPt );
  if( e3DList != x3Lst_tErr_NoErr )
    DebugPrint( ( "x3Lst error %d in NewCtrlPt.\n", e3DList ) );
  
  /* write it to the control point file. */
  if( ibWriteToFile )
    WriteVoxelToControlFile( ctrlPt );
}

void DeleteControlPoint ( xVoxelRef ipCtrlPt ) {
  
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  ctrlPt  = NULL;
  
  ctrlPt = ipCtrlPt;
  
  /* remove the item */
  e3DList = x3Lst_RemoveItem( gControlPointList, ctrlPt, (void**)&ctrlPt );
  if( e3DList != x3Lst_tErr_NoErr )
    goto error;
  
  /* delete it */
  xVoxl_Delete( &ctrlPt );
  
  goto cleanup;
  
 error:
  
  if( x3Lst_tErr_NoErr != e3DList
      && x3Lst_tErr_ItemNotInSpace != e3DList ) {
    DebugPrint( ( "x3Lst error %d in DeleteCtrlPt.\n", e3DList ) );
  }
  
 cleanup:
  return;
}

// ===========================================================================

// ========================================================= SELECTING REGIONS

tkm_tErr InitSelectionModule () {
  
  tkm_tErr   eResult = tkm_tErr_NoErr;
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  
  DebugEnterFunction( ("InitSelectionModule()") );
  
  DebugNote( ("Creating selction space") );
  e3DList = x3Lst_New( &gSelectedVoxels, NUM_UNDOS );
  DebugAssertThrowX( (x3Lst_tErr_NoErr == e3DList),
         eResult, tkm_tErr_Unrecoverable );
  
  DebugNote( ("Setting comparator in selection space") );
  e3DList = x3Lst_SetComparator( gSelectedVoxels, CompareVoxels );
  DebugAssertThrowX( (x3Lst_tErr_NoErr == e3DList),
         eResult, tkm_tErr_Unrecoverable );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void DeleteSelectionModule () {
  
  if ( NULL == gSelectedVoxels )
    return;
  
  DebugNote( ( "Deleting selection module.\n" ) );
  x3Lst_Delete( &gSelectedVoxels );
  
}

void AddVoxelsToSelection ( xVoxelRef iaAnaIdx, int inCount ) {
  
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  int        nVoxel  = 0;
  xVoxelRef  anaIdx  = NULL;
  
  if ( NULL == gSelectedVoxels )
    return;
  
  /* For each voxel we got... */
  for( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {
    
    /* Allocate a copy of the voxel */
    xVoxl_New( &anaIdx );
    xVoxl_Copy( anaIdx, &(iaAnaIdx[nVoxel]) );
    
    /* add the voxel to the ctrl pt space */
    e3DList = x3Lst_AddItem( gSelectedVoxels, anaIdx, anaIdx );
    if( e3DList != x3Lst_tErr_NoErr )
      DebugPrint( ( "x3Lst error %d in AddVoxelToSelection.\n", e3DList ) );
  }
}

void RemoveVoxelsFromSelection ( xVoxelRef iaAnaIdx, int inCount ) {
  
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  voxel   = NULL;
  int        nVoxel  = 0;
  xVoxel     where;
  
  if ( NULL == gSelectedVoxels )
    return;
  
  for( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {
    
    /* Copy the voxel into the location voxel. */
    xVoxl_Copy( &where, &(iaAnaIdx[nVoxel]) );
    
    voxel = &(iaAnaIdx[nVoxel]);
    
    /* Remove the item. we'll get back a ptr to the actual voxel we
       added. We may get an error if the voxel is out of bounds. If
       so, just continue. If we got a different error, print an error
       message.  */
    e3DList = x3Lst_RemoveItem( gSelectedVoxels, &where, (void**)&voxel );
    if( e3DList != x3Lst_tErr_NoErr ) {
      if( e3DList != x3Lst_tErr_ItemNotInSpace ) {
	DebugPrint( ( "x3Lst error %d in RemoveVoxelFromSelection. "
		      "voxel %d of %d = %d, %d, %d\n",
		      e3DList, nVoxel, inCount, xVoxl_ExpandInt( voxel ) ) );
      }
      continue;
    }
    
    /* delete it */
    xVoxl_Delete( &voxel );
  }
}


void ClearSelection () {
  
  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  
  if ( NULL == gSelectedVoxels )
    return;
  
  e3DList = x3Lst_Clear( gSelectedVoxels );
  if ( e3DList != x3Lst_tErr_NoErr )
    DebugPrint( ( "x3Lst error %d in ClearSelection\n",
      e3DList ) );
  
  MWin_RedrawAll( gMeditWindow );
}

void SaveSelectionToLabelFile ( char * isFileName ) {
  
  tkm_tErr  eResult       = tkm_tErr_NoErr;
  Volm_tErr  eVolume       = Volm_tErr_NoErr;
  char    sFileName[tkm_knPathLen] = "";
  LABEL*  pLabel       = NULL;
  LABEL_VERTEX* pVertex       = NULL;
  xVoxelRef  idx       = NULL;
  xVoxel  ras;
  xListRef  list       = NULL;
  x3Lst_tErr  e3DList       = x3Lst_tErr_NoErr;
  xList_tErr  eList       = xList_tErr_NoErr;
  int    nNumVoxels     = 0;
  int    nListCount     = 0;
  int    nVoxel       = 0;
  int    nPlane       = 0;
  int    eLabel       = 0;
  
  DebugEnterFunction( ("SaveSelectionToLabelFile ( isFileName=%s )",
           isFileName) );
  
  /* get the number of selected voxels we have, or the sum of the voxels
     in a particular plane. */
  for ( nPlane = 0; nPlane < gnAnatomicalDimensionZ; nPlane++ ) {
    
    DebugNote( ("Getting x items for plane %d", nPlane) );
    e3DList = x3Lst_GetItemsInXPlane( gSelectedVoxels, nPlane, &list );
    DebugAssertThrowX( (x3Lst_tErr_NoErr == e3DList),
           eResult, tkm_tErr_ErrorAccessingList );
    
    DebugNote( ("Getting number of items in list in plane %d", nPlane) );
    eList = xList_GetCount( list, &nListCount );
    DebugAssertThrowX( (xList_tErr_NoErr == eList),
           eResult, tkm_tErr_ErrorAccessingList );
    
    nNumVoxels += nListCount;
  }
  
  /* make the file name */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_Label, 
    sFileName, sizeof(sFileName) );
  
  /* allocate a label file with that number of voxels, our subject name,
     and the passed in label file name. */
  DebugNote( ("Allocating label with %d voxels") );  
  pLabel = LabelAlloc( nNumVoxels, NULL, sFileName );
  DebugAssertThrowX( (NULL != pLabel), eResult, tkm_tErr_CouldntAllocate );
  
  /* set the number of points in the label */
  pLabel->n_points = nNumVoxels;
  
  /* for each plane */
  nVoxel = 0;
  for ( nPlane = 0; nPlane < gnAnatomicalDimensionZ; nPlane++ ) {
    
    /* get the list for this x value. */
    DebugNote( ("Getting x items for plane %d", nPlane) );
    e3DList = x3Lst_GetItemsInXPlane( gSelectedVoxels, nPlane, &list );
    DebugAssertThrowX( (x3Lst_tErr_NoErr == e3DList),
           eResult, tkm_tErr_ErrorAccessingList );
    
    /* traverse the list */
    DebugNote( ("Resetting list position") );
    eList = xList_ResetPosition( list );
    while( (eList = xList_NextFromPos( list, (void**)&idx )) 
     != xList_tErr_EndOfList ) {
      
      if( idx ) {
  
  /* convert voxel to ras */
  eVolume = Volm_ConvertIdxToRAS(gAnatomicalVolume[tkm_tVolumeType_Main],
               idx, &ras );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* get a ptr the vertex in the label file. */
  DebugNote( ("Geting ptr to vertex %d", nVoxel) );
  pVertex = &(pLabel->lv[nVoxel]);
  
  /* set the vertex */
  DebugNote( ("Setting values of vertex %d", nVoxel) );
  pVertex->x = xVoxl_GetFloatX( &ras );
  pVertex->y = xVoxl_GetFloatY( &ras );
  pVertex->z = xVoxl_GetFloatZ( &ras );
  
  /* set the vno to -1, which is significant somewhere outside the realm
     of tkmedit. set stat value to something decent and deleted to not */
  pVertex->vno = -1;
  pVertex->stat = 0;
  pVertex->deleted = FALSE;
  
  /* inc our global count. */
  nVoxel++;
      }
    }
    
    DebugAssertThrowX( (eList == xList_tErr_EndOfList),
           eResult, tkm_tErr_ErrorAccessingList );
  }
  
  /* write the file */
  DebugNote( ("Writing label file to %s", sFileName) );
  eLabel = LabelWrite( pLabel, sFileName );
  DebugAssertThrowX( (NO_ERROR == eLabel),
         eResult, tkm_tErr_CouldntWriteFile );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != pLabel ) {
    DebugNote( ("Deleting label") );
    LabelFree( &pLabel );
  }
  
  DebugExitFunction;
}

tkm_tErr LoadSelectionFromLabelFile ( char* isFileName ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr  eVolume         = Volm_tErr_NoErr;
  char    sFileName[tkm_knPathLen]   = "";
  LABEL*  pLabel         = NULL;
  LABEL_VERTEX* pVertex         = NULL;
  int    nNumVoxels       = 0;
  int    nVoxel         = 0;
  xVoxel  ras;
  xVoxel  idx;
  char    sError[tkm_knErrStringLen] = "";
  
  DebugEnterFunction( ("LoadSelectionFromLabelFile ( isFileName=%s)",
           isFileName) );
  
  /* make the file name */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_Label, 
    sFileName, sizeof(sFileName) );
  
  /* read in the label */
  DebugNote( ("Reading label with LabelRead(%s)", sFileName));
  pLabel = LabelRead( NULL, sFileName );
  DebugAssertThrowX( (pLabel != NULL), eResult, tkm_tErr_CouldntLoadLabel );
  
  /* for each vertex in there... */
  nNumVoxels = pLabel->max_points;
  for( nVoxel = 0; nVoxel < nNumVoxels; nVoxel++ ) {
    
    /* get the vertex. */
    pVertex = &(pLabel->lv[nVoxel]);
    
    /* only process verticies that arn't deleted. */
    if ( !(pVertex->deleted) ) {
      
      /* transform from ras to voxel */
      xVoxl_SetFloat( &ras, pVertex->x, pVertex->y, pVertex->z );
      eVolume = Volm_ConvertRASToIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
              &ras, &idx );
      DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
       eResult, tkm_tErr_ErrorAccessingVolume );
      
      /* add to selection */
      AddVoxelsToSelection( &idx, 1 );
    }
  }    
  
  /* set the window's selection again to force a redraw. */
  MWin_SetSelectionSpace( gMeditWindow, -1, gSelectedVoxels );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading label %s", isFileName );
  tkm_DisplayError( sError,
		    tkm_GetErrorString(eResult),
		    "Tkmedit couldn't read the label you "
		    "specified. This could be because the format "
		    "wasn't valid or the file wasn't found." );
  EndDebugCatch;
  
  if( NULL != pLabel ) {
    DebugNote( ("Deleting label") );
    LabelFree( &pLabel );
  }
  
  DebugExitFunction;
  
  return eResult;
}

void GraphSelectedRegion () {
  
  xVoxelRef  voxel   = NULL;
  xListRef   list   = NULL;
  FunV_tErr  eFunctional = FunV_tErr_NoError;
  x3Lst_tErr e3DList   = x3Lst_tErr_NoErr;
  xList_tErr eList   = xList_tErr_NoErr;
  int       nPlane   = 0;
  
  // clear the functional display list.
  eFunctional = FunV_BeginSelectionRange( gFunctionalVolume );
  if( FunV_tErr_NoError != eFunctional )
    goto error;
  
  // get the ctrl pts in the list...
  for ( nPlane = 0; nPlane < gnAnatomicalDimensionZ; nPlane++ ) {
    
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
    DebugPrint( ( "xList error %d in GraphSelectedRegion.\n", eList ) );
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }
  
  if( e3DList != x3Lst_tErr_NoErr ) {
    DebugPrint( ( "x3Lst error %d in GraphSelectedRegion.\n", eList ) );
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }
  
  if( eFunctional != FunV_tErr_NoError ) {
    DebugPrint( ( "FunV error %d in GraphSelectedRegion.\n", 
      eFunctional ) );
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }
  
  
 cleanup:
  return;
}

void SelectVoxelsByFuncValue ( FunV_tFindStatsComp iCompare ) {
  
  xVoxel cursor;
  
  DebugEnterFunction( ("SelectVoxelsByFuncValue") );
  
  if( NULL == gMeditWindow ||
      NULL == gFunctionalVolume )
    return;
  
  /* get cursor */
  DebugNote( ("Getting cursor") );
  MWin_GetCursor( gMeditWindow, &cursor );
  
  /* select voxels */
  DebugNote( ("Selecting voxels around %d,%d,%d", xVoxl_ExpandInt( &cursor )));
  FunV_SelectAnaVoxelsByFuncValue( gFunctionalVolume, &cursor, iCompare );
  
  DebugExitFunction;
}


// ===========================================================================

// ============================================================= FILENAME MGMT

tkm_tErr SetSubjectHomeDirFromEnv ( char* isSubject ) {
  
  tkm_tErr eResult        = tkm_tErr_NoErr;
  char*     sEnvVar        = NULL;
  
  DebugEnterFunction( ("SetSubjectHomeDirFromEnv( isSubject=%s )",
           isSubject) );
  
  if (NULL != gsCommandLineSubjectsDir)
    sEnvVar = gsCommandLineSubjectsDir ;
  else
    sEnvVar = getenv( "SUBJECTS_DIR" );
  DebugAssertThrowX( (NULL != sEnvVar), 
         eResult, tkm_tErr_SUBJECTSDIRNotDefined );
  
  /* make the dir */
  xUtil_snprintf( gsSubjectHomeDir, sizeof(gsSubjectHomeDir),
      "%s/%s", sEnvVar, isSubject );
  
  /* send tcl update */
  tkm_SendTclCommand( tkm_tTclCommand_UpdateSubjectDirectory, gsSubjectHomeDir );
  
  DebugPrint( ("Set subject home dir to %s\n", gsSubjectHomeDir) );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr SetSubjectHomeDir ( char* isHomeDir ) {
  
  DebugEnterFunction( ("SetSubjectHomeDir( isHomeDir=%s )",
           isHomeDir) );
  
  xUtil_strncpy( gsSubjectHomeDir, isHomeDir, sizeof(gsSubjectHomeDir) );
  DebugPrint( ("Set subject home dir to %s\n", gsSubjectHomeDir) );
  
  DebugExitFunction;
  
  return tkm_tErr_NoErr;
}

tkm_tErr FindUserHomeDir () {
  
  tkm_tErr eResult        = tkm_tErr_NoErr;
  char*     psCurrentWorkingDir = NULL;
  
  DebugEnterFunction( ("FindUserHomeDir()") );
  
#ifdef Linux
  DebugNote( ("Finding pwd with getcwd") );
  psCurrentWorkingDir = getcwd(NULL,0);
#else
  DebugNote( ("Finding pwd with getenv(\"PWD\")") );
  psCurrentWorkingDir = getenv("PWD");
#endif
  
  DebugAssertThrowX( (NULL != psCurrentWorkingDir), 
         eResult, tkm_tErr_CouldntGetPWD );
  
  /* save cwd */
  xUtil_strncpy( gsUserHomeDir, psCurrentWorkingDir, sizeof(gsUserHomeDir) );
  
  DebugPrint( ( "Set user home dir to %s\n", gsUserHomeDir ) );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void MakeFileName ( char*     isInput,
		    tkm_tFileName  iType, 
		    char*         osCompleteFileName,
		    int         inDestSize ) {
  
  char sFileName[tkm_knPathLen] = "";
  
  DebugEnterFunction( ("MakeFileName( isInput=%s, iType=%d,"
           "osCompleteFileName=%s, inDestSize=%d )", 
           isInput, (int)iType, osCompleteFileName, inDestSize) );
  
  /* if null, return nothing */
  if( NULL == isInput )
    return;
  
  /* look at the first character */
  switch( isInput[0] ) {
    
    /* tilde, attach the subject home dir. */
  case '~':
    
    /* look at the second char. if it's a / or null use this subject. */
    if( isInput[1] == '/' ||
  isInput[1] == '\0' ) {
      
      xUtil_snprintf( sFileName, sizeof(sFileName),
          "%s%s", gsSubjectHomeDir, &(isInput[1]) );
      
      /* if it's not, treat the part after it as the subject name */
    } else {
      xUtil_snprintf( sFileName, sizeof(sFileName), 
          "%s/../%s", gsSubjectHomeDir, &(isInput[1]) );
    }
    
    break;
    
    /* period, use cwd nd then rest of input */
  case '.':
    xUtil_snprintf( sFileName, sizeof(sFileName),
        "%s%s", gsUserHomeDir, &(isInput[1]) );
    break;
    
    /* slash, it's a full path. */
  case '/':
    xUtil_strncpy( sFileName, isInput, sizeof(sFileName) );
    break;
    
    /* else, prepend subject home, then sub dir, then rest of file name. */
  default:
    if( gEnableFileNameGuessing ) {
      xUtil_snprintf( sFileName, sizeof(sFileName),
		      "%s/%s/%s", gsSubjectHomeDir, 
		      ksaFileNameSubDirs[iType], isInput );
    } else {
      xUtil_snprintf(sFileName,sizeof(sFileName), "./%s", isInput);
      /* xUtil_strncpy( sFileName, isInput, sizeof(sFileName) ); */
    }     
    break;
  }
  
  xUtil_strncpy( osCompleteFileName, sFileName, inDestSize );
  
  DebugExitFunction;
}

void ExtractSubjectName ( char* isDataSource,
        char* osSubjectName ) {
  
  int  nChar     = 0;
  int  nWordChar  = 0;
  char* sWord     = NULL;
  char  sName[STRLEN] = "";
  int  nLastSlash = 0;
  
  /* look for 'subjects' in the title */
  sWord = strstr( isDataSource, "/subjects/" );
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
    
    /* look for /mri in the title */
  } else if ( (sWord = strstr( isDataSource, "/mri" )) != NULL ){ 
    
    /* we're at the slash now. go to the slash before this one. while
       we're not at the beginning.. */
    while( sWord != isDataSource ) {
      sWord -= sizeof( char );
      if( *sWord == '/' )
  break;
    }
    
    /* inc past the slash and use the next part as the name */
    sWord += sizeof( char );
    nWordChar = 0;
    while( *sWord != '/' &&
     *sWord != '\0' ) {
      sName[nWordChar] = *sWord;
      nWordChar++;
      sWord += sizeof( char );
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
  
  int nChar   = 0;
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


/* ===================================================== General utilities */

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
      DebugPrint( ( "Cmd: %s\n", inCommand ) );
      DebugPrint( ( "\tCommand did not return OK.\n" ) );
      DebugPrint( ( "\tResult: %s\n", theInterp->result ) );
    }
    
  } else {
    
    /* cache the command so we can send it later */
    if( gNumCachedCommands < tkm_knMaxNumCachedCommands ) {
      strcpy( gCachedTclCommands[gNumCachedCommands++], inCommand );
    } else {
      DebugPrint( ( "Couldn't cache %s\n", inCommand ) );
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


tkm_tErr LoadVolume ( tkm_tVolumeType iType,
		      char*      isName ) {
  
  tkm_tErr eResult        = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  char     sPath[tkm_knPathLen]      = "";
  char*     pEnd          = NULL;
  char     sError[tkm_knErrStringLen]    = "";
  mriVolumeRef newVolume = NULL;
  
  DebugEnterFunction( ("LoadVolume( iType=%d,  isName=%s )", 
           (int)iType, isName) );
  
  DebugAssertThrowX( (NULL != isName), eResult, tkm_tErr_InvalidParameter );
  
  /* make a filename */
  DebugNote( ("Making filename from %s", isName) );
  MakeFileName( isName, tkm_tFileName_Volume, sPath, sizeof(sPath) );
  
  /* if there is a /COR- at the end, take it off. */
  DebugNote( ("Looking for COR- at end of file name") );
  pEnd = strstr( sPath, "/COR-" );
  if( NULL != pEnd )
    *pEnd = '\0';
  
  /* load the volume */
  DebugNote( ("Creating volume") );
  eVolume = Volm_New( &newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_CouldntReadVolume );
  
  DebugNote( ("Reading data into volume") );
  eVolume = Volm_ImportData( newVolume, sPath );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_CouldntReadVolume );
  
  /* if the volume exists, delete it */
  if( NULL != gAnatomicalVolume[iType] ) {
    UnloadDisplayTransform( iType );
    Volm_Delete( &gAnatomicalVolume[iType] );
  }
  
  /* save the new volume */
  gAnatomicalVolume[iType] = newVolume;

  /* show the tal coords and hide the ras coords */
  if (NULL != gAnatomicalVolume[iType]->mpMriValues->linear_transform)
    {
      DebugNote( ("Showing Tal coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowTalCoords, "1" );
      DebugNote( ("Hiding RAS coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowRASCoords, "0" );
    }
  else
    {
      DebugNote( ("Hiding Tal coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowTalCoords, "0" );
      DebugNote( ("Showing RAS coords") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowRASCoords, "1" );
    }
  
  /* save the volume size */
  DebugNote( ("Getting dimension of volume") );
  Volm_GetDimensions( gAnatomicalVolume[iType], &gnAnatomicalDimensionX, 
          &gnAnatomicalDimensionY, &gnAnatomicalDimensionZ  );
  
  
#if 1
  /* is this right (BRF)? Are they supposed to be screen dimensions? */
  gnAnatomicalDimensionX = gnAnatomicalDimensionY = gnAnatomicalDimensionZ = 256 ;
#else
  gnAnatomicalDimensionX = gAnatomicalVolume[iType]->mpMriValues->width ;
  gnAnatomicalDimensionY = gAnatomicalVolume[iType]->mpMriValues->height ;
  gnAnatomicalDimensionZ = gAnatomicalVolume[iType]->mpMriValues->depth ;
  {
    MATRIX *m_vox_to_ras ;
    VECTOR *v_vox, *v_ras ;
    int     r, c ;
    
    printf("w = %d, h = %d, d = %d\n",
     gAnatomicalVolume[iType]->mpMriValues->width,
     gAnatomicalVolume[iType]->mpMriValues->height,
     gAnatomicalVolume[iType]->mpMriValues->depth) ;
    
    m_vox_to_ras = MatrixInverse(gAnatomicalVolume[iType]->m_resample, NULL) ;
    for (r = 1 ; r <= 3 ; r++)
      {
  *MATRIX_RELT(m_vox_to_ras, r, 4) = 0 ;
  for (c = 1 ; c <= 3 ; c++)
    if (!FZERO(*MATRIX_RELT(m_vox_to_ras, r, c)))
      *MATRIX_RELT(m_vox_to_ras, r, c) = 1.0 ;
      }
    
    v_vox = VectorAlloc(4, 1) ;
    v_vox->rptr[4][1] = 1.0 ;
    V3_X(v_vox) = gAnatomicalVolume[iType]->mpMriValues->width ;
    V3_Y(v_vox) = gAnatomicalVolume[iType]->mpMriValues->height ;
    V3_Z(v_vox) = gAnatomicalVolume[iType]->mpMriValues->depth ;
    v_ras = MatrixMultiply(m_vox_to_ras, v_vox,NULL) ;
    gnAnatomicalDimensionX = fabs(V3_X(v_ras)) ;
    gnAnatomicalDimensionY = fabs(V3_Y(v_ras)) ;
    gnAnatomicalDimensionZ = fabs(V3_Z(v_ras)) ;
    MatrixFree(&m_vox_to_ras) ; MatrixFree(&v_vox) ; MatrixFree(&v_ras) ;
  }
  
#endif
  
  /*  if (Gdiag & DIAG_SHOW)
  printf("setting anatomical dimensions to %d, %d, %d\n",
   gnAnatomicalDimensionX, gnAnatomicalDimensionY, gnAnatomicalDimensionZ  );
  */

  /* set the default color scale */
  SetVolumeBrightnessAndContrast( iType, 
				  Volm_kfDefaultBrightness,
				  Volm_kfDefaultContrast );
  
  /* volume is clean */
  gbAnatomicalVolumeDirty[iType] = FALSE;
  
  /* set data in window */
  if( NULL != gMeditWindow ) {
    DebugNote( ("Setting volume in main window") );
    switch( iType ) {
    case tkm_tVolumeType_Main:  /* idx to ras xform always comes from main volume */
      MWin_SetVolume( gMeditWindow, -1, gAnatomicalVolume[iType],
          gnAnatomicalDimensionX, gnAnatomicalDimensionY, gnAnatomicalDimensionZ  );
      
      /* get a ptr to the idx to ras transform */
      DebugNote( ("Getting a pointer to the idx to RAS transform") );
      Volm_GetIdxToRASTransform( gAnatomicalVolume[iType],
         &gIdxToRASTransform );
      break;
    case tkm_tVolumeType_Aux:
      MWin_SetAuxVolume( gMeditWindow, -1, gAnatomicalVolume[iType],
       gnAnatomicalDimensionX, gnAnatomicalDimensionY, gnAnatomicalDimensionZ );
      break;
    default:
      break;
    }
  }
  

  gm_screen2ras = extract_i_to_r( gAnatomicalVolume[iType]->mpMriValues );
  gm_ras2screen = extract_r_to_i( gAnatomicalVolume[iType]->mpMriValues );


  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading volume %s", isName );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the volume you specified. "
        "This could be because the image format wasn't "
        "recognized, or it couldn't find the proper header, "
        "or the file(s) were unreadable." );
  
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


tkm_tErr UnloadVolume ( tkm_tVolumeType iType ) {

  tkm_tErr  eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("UnloadVolume( iType=%d )", (int)iType) );
  
  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumVolumeTypes),
		     eResult, tkm_tErr_InvalidParameter );

  /* Can't unload the main volume. */
  DebugAssertThrowX( (tkm_tVolumeType_Main != iType),
		     eResult, tkm_tErr_InvalidParameter );


  switch( iType ) {
  case tkm_tVolumeType_Aux:
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeOptions, "0" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeDirtyOptions, "0" );
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxTransformLoadedOptions, "0" );
    MWin_SetAuxVolume( gMeditWindow, -1, NULL, 
		       gnAnatomicalDimensionX, gnAnatomicalDimensionY,
		       gnAnatomicalDimensionZ );
    break;
  default:
    break;
  }

  /* If the volume exists, delete it */
  if( NULL != gAnatomicalVolume[iType] ) {
    UnloadDisplayTransform( iType );
    Volm_Delete( &gAnatomicalVolume[iType] );
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr LoadDisplayTransform ( tkm_tVolumeType iVolume,
        char*    isFileName ) {
  
  tkm_tErr  eResult           = tkm_tErr_NoErr;
  Volm_tErr eVolume           = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen]   = "";
  char      sError[tkm_knErrStringLen] = "";
  
  DebugEnterFunction( ("LoadDisplayTransform( iVolume=%d, isFileName=%s )",
           (int)iVolume, isFileName) );
  
  /* make a file name */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_VolumeTransform, 
    sFileName, sizeof(sFileName) );
  
  /* load the transform into our volume if we have it. */
  if( NULL != gAnatomicalVolume[iVolume] ) {
    eVolume = Volm_LoadDisplayTransform( gAnatomicalVolume[iVolume],
           sFileName );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume ),
           eResult, tkm_tErr_ErrorAccessingVolume );
  }
  
  /* enable the menu options */
  switch( iVolume ) {
  case tkm_tVolumeType_Main:
    tkm_SendTclCommand( tkm_tTclCommand_ShowMainTransformLoadedOptions, "1" );
    break;
  case tkm_tVolumeType_Aux:
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxTransformLoadedOptions, "1" );
    break;
  default:
    break;
  }
  
  /* big redraw */  
  if( NULL != gMeditWindow ) 
    MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading transform %s", isFileName );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the transform you specified. "
        "This could be because the file format wasn't "
        "recognized, or the file was unreadable." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr UnloadDisplayTransform ( tkm_tVolumeType iVolume ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("UnloadDisplayTransform( iVolume=%d )", (int)iVolume) );
  
  /* unload the transform in the volume if we have it */
  if( NULL != gAnatomicalVolume[iVolume] ) {
    eVolume = Volm_UnloadDisplayTransform( gAnatomicalVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume ),
           eResult, tkm_tErr_ErrorAccessingVolume );
  }
  
  /* disable the menu options */
  switch( iVolume ) {
  case tkm_tVolumeType_Main:
    tkm_SendTclCommand( tkm_tTclCommand_ShowMainTransformLoadedOptions, "0" );
    break;
  case tkm_tVolumeType_Aux:
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxTransformLoadedOptions, "0" );
    break;
  default:
    break;
  }
  
  /* big redraw */  
  if( NULL != gMeditWindow ) 
    MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr SaveVolume ( tkm_tVolumeType iVolume,
		      char*           isPath ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  tBoolean  bDirty         = FALSE;
  char      sFileName[tkm_knPathLen] = "";
  char*      psFileName         = NULL;
  
  DebugEnterFunction( ("SaveVolume( iVolume=%d, isPath=%s )", 
		       iVolume, isPath) );
  
  if( NULL == gAnatomicalVolume[iVolume] ) 
    DebugGotoCleanup;
  
  /* if we have editing disabled, return */
  if( !editflag ) {
    tkm_DisplayError( "Saving Volume", 
		      "Couldn't save volume",
		      "This session was started in read-only mode. You cannot "
		      "save changes." );
    DebugGotoCleanup;
  }
  
  /* if we haven't been edited, return */
  eResult = IsVolumeDirty( iVolume, &bDirty );
  if( !bDirty ) {
    DebugPrint( ("SaveVolume called when volume is clean!\n") );
    DebugGotoCleanup;
  }
  
  /* make a file name */
  if( NULL != isPath ) {
    DebugNote( ("Making file name from %s", isPath) );
    MakeFileName( isPath, tkm_tFileName_Volume, sFileName, sizeof(sFileName) );
    psFileName = sFileName;
  } 
  
  /* write the main anatomical volume */
  eVolume = Volm_ExportNormToCOR( gAnatomicalVolume[iVolume], 
          psFileName );
  if( Volm_tErr_NoErr != eVolume ) {
    tkm_DisplayError( "Saving Volume", 
          "Couldn't save volume",
          "The files could not be written. Check to make sure the "
          "destination directory exists, that you have "
          "permission to write to it, and that there is "
          "sufficient disk space." );
  }
  
  /* set our edited flags null */
  eResult = SetVolumeDirty( iVolume, FALSE );
  DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
  
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void SnapshotVolume () {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("SnapshotVolume ()") );
  
  /* tell the volume to save a snapshot */
  DebugNote( ("Saving snapshot") );
  eVolume = Volm_SaveToSnapshot( gAnatomicalVolume[tkm_tVolumeType_Main] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void RestoreVolumeFromSnapshot () {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("RestoreVolumeFromSnapshot ()") );
  
  /* tell the volume to save a snapshot */
  DebugNote( ("Restoring snapshot") );
  eVolume = Volm_RestoreFromSnapshot(gAnatomicalVolume[tkm_tVolumeType_Main]);
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* big redraw */
  MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

tkm_tErr SetVolumeDirty ( tkm_tVolumeType iVolume, tBoolean ibDirty ) {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("SetVolumeDirty ( iVolume=%d, ibDirty=%d )",
           (int)iVolume, (int)ibDirty) );
  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
         eResult, tkm_tErr_InvalidParameter );
  
  /* set the dirty flag */
  gbAnatomicalVolumeDirty[iVolume] = ibDirty;
  
  /* updatex tcl */
  switch( iVolume ) {
  case tkm_tVolumeType_Main:
    tkm_SendTclCommand( tkm_tTclCommand_ShowVolumeDirtyOptions, 
			ibDirty?"1":"0" );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeDirty, ibDirty?"1":"0" );
    break;
  case tkm_tVolumeType_Aux:
    tkm_SendTclCommand( tkm_tTclCommand_ShowAuxVolumeDirtyOptions, 
			ibDirty?"1":"0" );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateAuxVolumeDirty,ibDirty?"1":"0" );
    break;
  default:
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr IsVolumeDirty ( tkm_tVolumeType iVolume, tBoolean* obDirty ) {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("IsVolumeDirty ( iVolume=%d, obDirty=%p )",
		       (int)iVolume, obDirty) );
  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
		     eResult, tkm_tErr_InvalidParameter );
  
  /* get the dirty flag */
  *obDirty = gbAnatomicalVolumeDirty[iVolume];
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void SetVolumeBrightnessAndContrast ( tkm_tVolumeType iVolume,
              float        ifBrightness,
              float        ifContrast ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("SetVolumeBrightnessAndContrast ( iVolume=%d, "
		       "ifBrightness=%.2f, ifContrast=%.2f )",
		       (int)iVolume, ifBrightness, ifContrast) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes), 
		     eResult, tkm_tErr_InvalidParameter );

  /* save these values. */
  gfaBrightness[iVolume] = ifBrightness;
  gfaContrast[iVolume] = ifContrast;

  /* set the b/c in the volume. */
  DebugNote( ("Setting brightness and contrast") );
  eVolume = Volm_SetBrightnessAndContrast( gAnatomicalVolume[iVolume], 
					   gfaBrightness[iVolume],
					   gfaContrast[iVolume] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* update the tcl window */
  SendBrightnessAndContrastUpdate( iVolume );
  
  /* big redraw */
  MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void SendBrightnessAndContrastUpdate ( tkm_tVolumeType iVolume ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  char      sTclArguments[tkm_knTclCmdLen] = "";
  
  DebugEnterFunction( ("SendBrightnessAndContrastUpdate ( iVolume=%d )",
		       (int)iVolume) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes), 
		     eResult, tkm_tErr_InvalidParameter );

  /* update the tcl window */
  DebugNote( ("Sending color scale update to tcl window") );
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %f %f",
		  (int)iVolume, gfaBrightness[iVolume], gfaContrast[iVolume] );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeColorScale, sTclArguments );
  
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void ThresholdVolume ( int    inLevel,
		       tBoolean   ibAbove,
		       int        inNewLevel ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("ThresholdVolume ( inLevel=%d, ibAbove=%d, "
           "inNewLevel=%d)", inLevel, (int)ibAbove, inNewLevel) );
  
  /* tell the volume to threshold */
  DebugNote( ("Thresholding volume") );
  eVolume = Volm_Threshold( gAnatomicalVolume[tkm_tVolumeType_Main],
          (Volm_tValue)inLevel, ibAbove, 
          (Volm_tValue)inNewLevel);
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* volume is dirty. */
  SetVolumeDirty( tkm_tVolumeType_Main, TRUE );
  
  /* big redraw */
  MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void FlipVolume ( mri_tOrientation  iAxis ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("FlipVolume ( iAxis=%d )", iAxis) );
  
  /* tell the volume to flip */
  DebugNote( ("Flipping volume") );
  eVolume = Volm_Flip( gAnatomicalVolume[tkm_tVolumeType_Main], iAxis );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* volume is dirty. */
  SetVolumeDirty( tkm_tVolumeType_Main, TRUE );
  
  /* big redraw */
  MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void RotateVolume ( mri_tOrientation  iAxis,
        float      ifDegrees ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("RotateVolume ( iAxis=%d, ifDegrees=%.2f )",
           iAxis, ifDegrees) );
  
  /* tell the volume to threshold */
  DebugNote( ("Rotating volume") );
  eVolume = Volm_Rotate( gAnatomicalVolume[tkm_tVolumeType_Main],
       iAxis, ifDegrees );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* volume is dirty. */
  SetVolumeDirty( tkm_tVolumeType_Main, TRUE );
  
  /* big redraw */
  MWin_RedrawAll( gMeditWindow );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

void EditAnatomicalVolumeInRangeArray ( tkm_tVolumeType iVolume,
					xVoxelRef       iaVoxel, 
					int             inCount,
					Volm_tValue     inLow,
					Volm_tValue     inHigh, 
					Volm_tValue     inNewValue ) {
  
  int          nVoxel   = 0;
  float        value    = 0;
  Volm_tValue  newValue = 0;

  if( NULL == gAnatomicalVolume[iVolume] ) {
    return;
  }
  
  /* for each voxel we got... */
  for( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {

    /* get the value at this point. */
    Volm_GetValueAtIdx( gAnatomicalVolume[iVolume], 
			&(iaVoxel[nVoxel]), &value );
    newValue = (int)value;
  
    /* if it's in the range... */
    if( value >= inLow && value <= inHigh ) {
      
      /* set a new value */
      newValue = inNewValue;
    }
    
    /* if values are different, set and add to undo list. */
    if( value != newValue ) {
      Volm_SetValueAtIdx( gAnatomicalVolume[iVolume], 
			  &(iaVoxel[nVoxel]), newValue );
      
      switch( iVolume ) {
      case tkm_tVolumeType_Main:
	/* if this is an edit on the main volume, add it to the undo list
	   using the EditAnatomicalVolume function and also add it to the
	   undo volume. */
	AddVoxelAndValueToUndoList( EditAnatomicalVolume, 
				    &(iaVoxel[nVoxel]), value ); 
	AddAnaIdxAndValueToUndoVolume( &(iaVoxel[nVoxel]), value ); 
	break;
	
      case tkm_tVolumeType_Aux:
	/* if this is an edit on the aux volume, add it to the undo list
	   using the EditAuxAnatomicalVolume function. */
	AddVoxelAndValueToUndoList( EditAuxAnatomicalVolume, 
				    &(iaVoxel[nVoxel]), value ); 
	break;

      default:
      }
    }
  }
  
  /* volume is dirty. */
  SetVolumeDirty( iVolume, TRUE );
}

void SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
				 int             iAnaX0,
				 int             iAnaX1,
				 int             iAnaY0,
				 int             iAnaY1,
				 int             iAnaZ0,
				 int             iAnaZ1,
				 float           iNewValue ) {
  
  tkm_tErr   eResult = tkm_tErr_NoErr;
  Volm_tErr  eVolume = Volm_tErr_NoErr;
  xVoxel     begin;
  xVoxel     end;
  xVoxel     cur;

  DebugEnterFunction( ("SetAnatomicalVolumeRegion( iVolume=%d, "
		       "iAnaX0=%d, iAnaY0=%d, iAnaZ0=%d, "
		       "iAnaX1=%d, iAnaY1=%d, iAnaZ1=%d )",
		       iVolume, iAnaX0, iAnaY0, iAnaZ0,
		       iAnaX1, iAnaY1, iAnaZ1) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
		     eResult, tkm_tErr_InvalidParameter );
  
  /* Make sure we got a good range. */
  xVoxl_Set( &begin, 
	     MIN( iAnaX0, iAnaX1 ),
	     MIN( iAnaY0, iAnaY1 ),
	     MIN( iAnaZ0, iAnaZ1 ) );
  eVolume = Volm_VerifyIdxInMRIBounds( gAnatomicalVolume[iVolume], &begin );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
		     eResult, tkm_tErr_InvalidParameter );
  xVoxl_Set( &end, 
	     MAX( iAnaX0, iAnaX1 ),
	     MAX( iAnaY0, iAnaY1 ),
	     MAX( iAnaZ0, iAnaZ1 ) );
  eVolume = Volm_VerifyIdxInMRIBounds( gAnatomicalVolume[iVolume], &end );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
		     eResult, tkm_tErr_InvalidParameter );

  /* step through the volume and set everything to this value. */
  DebugNote( ("Setting volume values") );
  xVoxl_Copy( &cur, &begin );
  while( xVoxl_IncrementWithMinsUntilLimits( &cur, 
					     xVoxl_GetX(&begin), 
					     xVoxl_GetY(&begin), 
					     xVoxl_GetX(&end), 
					     xVoxl_GetY(&end), 
					     xVoxl_GetZ(&end) ) ) {
    Volm_SetValueAtIdx_( gAnatomicalVolume[iVolume], &cur, iNewValue );
  }
  
  /* volume is dirty. */
  SetVolumeDirty( iVolume, TRUE );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
}

int EditAnatomicalVolume ( xVoxelRef iAnaIdx, int inValue ) {
  
  tkm_tErr    eResult = tkm_tErr_NoErr;
  Volm_tErr   eVolume = Volm_tErr_NoErr;
  float        value   = 0;
  
  DebugEnterFunction( ("EditAnatomicalVolume( iAnaIdx=%d,%d,%d, inValue=%d)",
		       xVoxl_ExpandInt(iAnaIdx), inValue) );
  
  /* get the current value so we can return it. set the new value */
  eVolume = Volm_GetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Main], 
				iAnaIdx, &value );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingVolume );
  eVolume = Volm_SetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Main], 
				iAnaIdx, inValue );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingVolume );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return (int)value;
}

int EditAuxAnatomicalVolume ( xVoxelRef iAnaIdx, int inValue ) {
  
  tkm_tErr    eResult = tkm_tErr_NoErr;
  Volm_tErr   eVolume = Volm_tErr_NoErr;
  float        value   = 0;
  
  DebugEnterFunction(("EditAuxAnatomicalVolume( iAnaIdx=%d,%d,%d, inValue=%d)",
		       xVoxl_ExpandInt(iAnaIdx), inValue) );
  
  /* get the current value so we can return it. set the new value */
  eVolume = Volm_GetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Aux], 
				iAnaIdx, &value );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingVolume );
  eVolume = Volm_SetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Aux], 
				iAnaIdx, inValue );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingVolume );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return (int)value;
}

void ConvertRASToAnaIdx ( xVoxelRef iRAS,
        xVoxelRef oAnaIdx ) {
  
  Volm_ConvertRASToIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
      iRAS, oAnaIdx );
}

void ConvertAnaIdxToRAS ( xVoxelRef iAnaIdx,
        xVoxelRef oRAS ) {
  
  Volm_ConvertIdxToRAS( gAnatomicalVolume[tkm_tVolumeType_Main],
      iAnaIdx, oRAS );
}
// ========================================================= FUNCTIONAL VOLUME

tkm_tErr LoadFunctionalOverlay( char* isPathAndStem, 
        char* isOffsetPathAndStem,
        char* isRegistration ) {
  
  tkm_tErr  eResult       = tkm_tErr_NoErr;
  char      sPathAndStem[tkm_knPathLen]   = "";
  char      sOffsetPathAndStem[tkm_knPathLen]  = "";
  char      sRegistration[tkm_knPathLen] = "";
  char*      psOffsetPathAndStem     = NULL;
  char*      psRegistration     = NULL;
  FunV_tErr eFunctional       = FunV_tErr_NoError;
  char      sError[tkm_knErrStringLen]   = "";
  
  DebugEnterFunction( ("LoadFunctionalOverlay( isPathAndStem=%s, "
           "isOffsetPathAndStem=%s, isRegistration=%s )", 
           isPathAndStem, isOffsetPathAndStem, isRegistration) );
  
  /* make our path and stem filename. if we have a registraltion file name,
     make that too. */
  DebugNote( ("Making file name from %s", isPathAndStem) );
  MakeFileName( isPathAndStem, tkm_tFileName_Functional, 
    sPathAndStem, sizeof(sPathAndStem) );
  if( NULL != isOffsetPathAndStem ) {
    DebugNote( ("Making file name from %s", isOffsetPathAndStem) );
    MakeFileName( isOffsetPathAndStem, tkm_tFileName_Functional, 
      sOffsetPathAndStem, sizeof(sOffsetPathAndStem) );
    psOffsetPathAndStem = sOffsetPathAndStem;
  } else {
    psOffsetPathAndStem = NULL;
  }
  if( NULL != isRegistration ) {
    DebugNote( ("Making file name from %s", isRegistration) );
    MakeFileName( isRegistration, tkm_tFileName_Functional, 
      sRegistration, sizeof(sRegistration) );
    psRegistration = sRegistration;
  } else {
    psRegistration = NULL;
  }
  
  /* attempt to load. */
  DebugNote( ("Loading overlay") );
  eFunctional = FunV_LoadOverlay( gFunctionalVolume, 
          sPathAndStem,
          psOffsetPathAndStem,
          psRegistration );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
         eResult, tkm_tErr_CouldntLoadOverlay );
  
  /* turn overlay display on */
  MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_FunctionalOverlay, 
           TRUE );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading functional overlay %s", 
      isPathAndStem );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the overlay volume you "
        "specified. This could be because the volume "
        "wasn't found or was unreadable, or because a "
        "valid header type couldn't be find, or a "
        "registration file couldn't be found or opened." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr LoadFunctionalTimeCourse( char* isPathAndStem,
           char* isOffsetPathAndStem,
           char* isRegistration ) {
  
  tkm_tErr  eResult            = tkm_tErr_NoErr;
  char      sPathAndStem[tkm_knPathLen]        = "";
  char      sOffsetPathAndStem[tkm_knPathLen] = "";
  char      sRegistration[tkm_knPathLen]      = "";
  char*      psOffsetPathAndStem          = NULL;
  char*      psRegistration          = NULL;
  FunV_tErr eFunctional            = FunV_tErr_NoError;
  char      sError[tkm_knErrStringLen]        = "";
  
  DebugEnterFunction( ("LoadFunctionalTimeCourse( isPathAndStem=%s, "
           "isOffsetPathAndStem=%s, isRegistration=%s )", 
           isPathAndStem, isOffsetPathAndStem, isRegistration) );
  
  /* make our path and stem filename. if we have a registration file name,
     make that too. */
  DebugNote( ("Making file name from %s", isPathAndStem) );
  MakeFileName( isPathAndStem, tkm_tFileName_Functional, 
    sPathAndStem, sizeof(sPathAndStem) );
  if( NULL != isOffsetPathAndStem ) {
    DebugNote( ("Making file name from %s", isOffsetPathAndStem) );
    MakeFileName( isOffsetPathAndStem, tkm_tFileName_Functional, 
      sOffsetPathAndStem, sizeof(sOffsetPathAndStem) );
    psOffsetPathAndStem = sOffsetPathAndStem;
  } else {
    psOffsetPathAndStem = NULL;
  }
  if( NULL != isRegistration ) {
    DebugNote( ("Making file name from %s", isRegistration) );
    MakeFileName( isRegistration, tkm_tFileName_Functional, 
      sRegistration, sizeof(sRegistration) );
    psRegistration = sRegistration;
  } else {
    psRegistration = NULL;
  }
  
  /* attempt to load. */
  DebugNote( ("Loading time course") );
  eFunctional = FunV_LoadTimeCourse( gFunctionalVolume, 
             sPathAndStem, 
             psOffsetPathAndStem,
             psRegistration );
  
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
         eResult, tkm_tErr_CouldntLoadOverlay );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading functional time course %s", 
      isPathAndStem );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the time course volume you "
        "specified. This could be because the volume "
        "wasn't found or was unreadable, or because a "
        "valid header type couldn't be find, or a "
        "registration file couldn't be found or opened." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr SmoothOverlayData ( float ifSigma ) {
  
  tkm_tErr  eResult            = tkm_tErr_NoErr;
  FunV_tErr eFunctional            = FunV_tErr_NoError;
  char      sError[tkm_knErrStringLen]        = "";
  
  DebugEnterFunction( ("SmoothOverlayData( ifSigma=%f )", ifSigma) );
  
  /* make sure we have data */
  DebugAssertThrowX( (NULL != gFunctionalVolume), 
         eResult, tkm_tErr_OverlayNotLoaded);
  
  /* smooth the data */
  printf("smoothing with sigma = %2.2f...\n", ifSigma) ;
  eFunctional = FunV_SmoothOverlayData( gFunctionalVolume, ifSigma );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
         eResult, tkm_tErr_ErrorAccessingFunctionalVolume );
  
  UpdateAndRedraw();
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Smoothing overlay data" );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't smooth the functional overlay data." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


void RotateOverlayRegistration ( float ifDegrees, 
         char  isAxis ) {
  
  mri_tOrientation orientation;
  xVoxel     cursor;
  
  MWin_GetOrientation ( gMeditWindow, &orientation );
  
  MWin_GetCursor( gMeditWindow, &cursor );
  
  switch( isAxis ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              ifDegrees, tAxis_X, &cursor );
      break;
    case mri_tOrientation_Horizontal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              -ifDegrees, tAxis_X, &cursor );
      break;
    case mri_tOrientation_Sagittal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              -ifDegrees, tAxis_Y, &cursor );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              ifDegrees, tAxis_Z, &cursor );
      break;
    case mri_tOrientation_Horizontal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              ifDegrees, tAxis_Y, &cursor );
      break;
    case mri_tOrientation_Sagittal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              ifDegrees, tAxis_Z, &cursor );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              -ifDegrees, tAxis_Y, &cursor );
      break;
    case mri_tOrientation_Horizontal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              ifDegrees, tAxis_Z, &cursor );
      break;
    case mri_tOrientation_Sagittal:
      FunV_RotateOverlayRegistration( gFunctionalVolume, 
              -ifDegrees, tAxis_X, &cursor );
      break;
    default: goto cleanup;
      break;
    }
    break;
  }
  
  UpdateAndRedraw();
  
 cleanup:
  return;
  
}

void TranslateOverlayRegisgtration ( float ifDistance, 
             char  isDirection ) {
  
  mri_tOrientation orientation;
  
  MWin_GetOrientation ( gMeditWindow, &orientation );
  
  switch( isDirection ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           -ifDistance, tAxis_Y );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           ifDistance, tAxis_Z );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           -ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Horizontal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           -ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      FunV_TranslateOverlayRegistration( gFunctionalVolume,
           -ifDistance, tAxis_X );
      break;
    default: goto cleanup;
      break;
    }
    break;
  }
  
  UpdateAndRedraw();
  
 cleanup:
  return;
  
}

void ScaleOverlayRegisgtration ( float ifFactor, 
         char  isDirection ) {
  
  mri_tOrientation orientation;
  
  MWin_GetOrientation ( gMeditWindow, &orientation );
  
  switch( isDirection ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Y );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Z );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Y );
      
    case mri_tOrientation_Horizontal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      FunV_ScaleOverlayRegistration( gFunctionalVolume,
             ifFactor, tAxis_X );
      break;
    default: goto cleanup;
      break;
    }
    break;
  default:
    goto cleanup;
  }
  
  UpdateAndRedraw();
  
 cleanup:
  return;
  
}



// ============================================================ SEGMENTATION 

tkm_tErr NewSegmentationVolume ( tkm_tSegType    iVolume,
				 tkm_tVolumeType iFromVolume,
				 char*           isColorFileName ) {
  
  tkm_tErr     eResult                      = tkm_tErr_NoErr;
  Volm_tErr    eVolume                      = Volm_tErr_NoErr;
  mriVolumeRef newVolume                    = NULL;
  mriVolumeRef newChangedVolume             = NULL;
  char         sError[tkm_knErrStringLen]   = "";

  DebugEnterFunction( ("LoadSegmentationVolume( iVolume=%d, iFromVolume=%d, "
		       "isColorFileName=%s )", (int)iVolume,
		       (int)iFromVolume, isColorFileName) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isColorFileName),
		     eResult, tkm_tErr_InvalidParameter );
  

  /* Create the new volume from the existing anatomical. Set it to all
     zeroes. */
  DebugNote( ("Creating new segmentation volume") );
  eVolume = Volm_New( &newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newVolume, gAnatomicalVolume[iFromVolume] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_SetAllValues( newVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Try to load the color table. */
  DebugNote( ("Loading color table.") );
  eResult = LoadSegmentationColorTable( iVolume, isColorFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  
  /* allocate flag volume from the new segmentation volume. set
     everything in it to zero */
  DebugNote( ("Creating segmentation flag volume") );
  eVolume = Volm_New( &newChangedVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newChangedVolume, newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_CouldntLoadSegmentation );
  eVolume = Volm_SetAllValues( newChangedVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  


  /* free existing segmentations if present. */
  if( NULL != gSegmentationVolume[iVolume] ) {
    DebugNote( ("Deleting existing roi group") );
    eVolume = Volm_Delete( &gSegmentationVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
    
    DebugNote( ("Deleting existing flag volume") );
    eVolume = Volm_Delete( &gSegmentationChangedVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
           eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }
  

  /* Save these volumes. */
  gSegmentationVolume[iVolume]        = newVolume;
  gSegmentationChangedVolume[iVolume] = newChangedVolume;

  /* enable our segmentation options and set segmentation volume in window */
  if( gMeditWindow ) {
    if( tkm_tSegType_Main == iVolume ) {
      DebugNote( ("Enabling segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowSegmentationOptions, "1" );
    } else {
      DebugNote( ("Enabling ayx segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowAuxSegmentationOptions, "1" );
    }
  
    DebugNote( ("Setting segmentation in main window") );
    MWin_SetSegmentationVolume( gMeditWindow, iVolume, 
				-1, gSegmentationVolume[iVolume] );
    MWin_SetSegmentationColorTable( gMeditWindow, iVolume, 
				    -1, gColorTable[iVolume] );
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr || tkm_tErr_CouldntLoadColorTable,
		   tkm_GetErrorString );
  
  if( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError), "Creating Segmentation" );
    tkm_DisplayError( sError, tkm_GetErrorString(eResult),
		      "Tkmedit couldn't create a segmentation." );
  }

  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr LoadSegmentationVolume ( tkm_tSegType iVolume,
				  char*        isVolumeDirWithPrefix,
				  char*        isColorFileName ) {
  
  tkm_tErr     eResult         = tkm_tErr_NoErr;
  Volm_tErr    eVolume         = Volm_tErr_NoErr;
  mriVolumeRef newVolume                    = NULL;
  mriVolumeRef newChangedVolume             = NULL;
  char         sSegmentationFileName[tkm_knPathLen] = "";
  char         sError[tkm_knErrStringLen]     = "";

  DebugEnterFunction( ("LoadSegmentationVolume( iVolume=%d, "
		       "isVolumeDirWithPrefix=%s, isColorFileName=%s )", 
		       iVolume, isVolumeDirWithPrefix, isColorFileName) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isVolumeDirWithPrefix),
         eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isColorFileName),
         eResult, tkm_tErr_InvalidParameter );
  
  /* make a file name */
  DebugNote( ("Making file name from %s", isVolumeDirWithPrefix) );
  MakeFileName( isVolumeDirWithPrefix, tkm_tFileName_Segmentation, 
    sSegmentationFileName, sizeof(sSegmentationFileName) );
    
  /* read in segmentation volume */
  DebugNote( ("Creating segmentation overlay") );
  eVolume = Volm_New( &newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  DebugNote( ("Importing segmentation into roi group") );
  eVolume = Volm_ImportData( newVolume, sSegmentationFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_CouldntLoadSegmentation );
  
  /* Try to load the color table. */
  DebugNote( ("Loading color table.") );
  eResult = LoadSegmentationColorTable( iVolume, isColorFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* allocate changed volume */
  DebugNote( ("Creating segmentation flag volume") );
  eVolume = Volm_New( &newChangedVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newChangedVolume, newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_CouldntLoadSegmentation );
  eVolume = Volm_SetAllValues( newChangedVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* free existing roi group if present. */
  if( NULL != gSegmentationVolume[iVolume] ) {
    DebugNote( ("Deleting existing roi group") );
    eVolume = Volm_Delete( &gSegmentationVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
           eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
    
    DebugNote( ("Deleting existing flag volume") );
    eVolume = Volm_Delete( &gSegmentationChangedVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }
  
  /* Save these volumes. */
  gSegmentationVolume[iVolume]        = newVolume;
  gSegmentationChangedVolume[iVolume] = newChangedVolume;

  /* enable our segmentation options and set segmentation volume in window */
  if( gMeditWindow ) {
    if( tkm_tSegType_Main == iVolume ) {
      DebugNote( ("Enabling segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowSegmentationOptions, "1" );
    } else {
      DebugNote( ("Enabling ayx segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowAuxSegmentationOptions, "1" );
    }

    DebugNote( ("Setting segmentation in main window") );
    MWin_SetSegmentationVolume( gMeditWindow, iVolume, 
		      -1, gSegmentationVolume[iVolume] );
    MWin_SetSegmentationColorTable( gMeditWindow, iVolume, 
				    -1, gColorTable[iVolume] );
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr || tkm_tErr_CouldntLoadColorTable,
		   tkm_GetErrorString );
  
  if( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError),
		    "Loading Segmentation %s", isVolumeDirWithPrefix );
    tkm_DisplayError( sError,
		      tkm_GetErrorString(eResult),
		      "Tkmedit couldn't read the segmentation you "
		      "specified. This could be because the segmentation "
		      "volume wasn't a valid COR volume directory." );
  }

  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void SaveSegmentationVolume ( tkm_tSegType iVolume,
			      char*        isVolumeDirWithPrefix ) {
  
  tkm_tErr  eResult         = tkm_tErr_NoErr;
  char      sError[tkm_knErrStringLen] = "";
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sSegmentationFileName[tkm_knPathLen] = "";
  char*     psSegmentationFileName     = NULL;
  
  DebugEnterFunction( ("SaveSegmentationVolume ( iVolume=%d, "
		       "isVolumeDirWithPrefix=%s )", (int)iVolume,
		       isVolumeDirWithPrefix) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume]),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* make a file name if they gave us one. */
  if( NULL != isVolumeDirWithPrefix ) {
    MakeFileName( isVolumeDirWithPrefix, tkm_tFileName_Segmentation, 
		  sSegmentationFileName, sizeof(sSegmentationFileName) );
    psSegmentationFileName = sSegmentationFileName;
  } else {
    psSegmentationFileName = NULL;
  }
  
  /* Export the thing to a COR volume. */
  OutputPrint "Saving... " EndOutputPrint;
  DebugNote( ("Exporting to COR") );
  eVolume = Volm_ExportNormToCOR( gSegmentationVolume[iVolume], 
				  psSegmentationFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
		     eResult, tkm_tErr_CouldntWriteFile );
  OutputPrint "done.\n" EndOutputPrint;
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  if( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError),
		    "Saving Segmentation %s", isVolumeDirWithPrefix );
    tkm_DisplayError( sError,
		      tkm_GetErrorString(eResult),
		      "Tkmedit couldn't save the segmentation you "
		      "specified. This could be because the segmentation "
		      "volume wasn't loaded or the destination directory "
		      "wasn't valid." );
  }

  EndDebugCatch;
  
  DebugExitFunction;
}

tkm_tErr ExportChangedSegmentationVolume ( tkm_tSegType iVolume,
					   char*        isVolumeDir ) {
  
  tkm_tErr                      eResult                    = tkm_tErr_NoErr;
  Volm_tErr                     eVolume                    = Volm_tErr_NoErr;
  mriVolumeRef                  newVolume                  = NULL;
  char                          sSegmentationFileName[tkm_knPathLen] = "";
  tkm_tExportSegmentationParams params;
  
  DebugEnterFunction( ("ExportChangedSegmentationVolume ( iVolume=%d, "
		       "isVolumeDir=%s )", (int)iVolume, isVolumeDir) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isVolumeDir),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume]),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  DebugAssertThrowX( (NULL != gSegmentationChangedVolume[iVolume]),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* make a file name. */
  if( NULL != isVolumeDir ) {
    MakeFileName( isVolumeDir, tkm_tFileName_Segmentation, 
		  sSegmentationFileName, sizeof(sSegmentationFileName) );
  }
  
  OutputPrint "Saving... " EndOutputPrint;

  /* Create the new volume from the existing segmentation. Set it to all
     zeroes. */
  DebugNote( ("Creating new segmentation volume") );
  eVolume = Volm_New( &newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newVolume, gSegmentationVolume[iVolume] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_SetAllValues( newVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Fill out the params. */
  params.mSrcVolume = gSegmentationVolume[iVolume];
  params.mDestVolume = newVolume;

  /* Export the thing to a COR volume. */
  DebugNote( ("Running compare") );
  eVolume = 
    Volm_VisitAllVoxels( gSegmentationChangedVolume[iVolume],
			 &SetChangedSegmentationValue, (void*)&params );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Now save it. */
  eVolume = Volm_ExportNormToCOR( newVolume, sSegmentationFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), 
		     eResult, tkm_tErr_CouldntWriteFile );

  /* And get rid of the volume. */
  Volm_Delete( &newVolume );

  OutputPrint "done.\n" EndOutputPrint;
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr ImportSurfaceAnnotationToSegmentation ( tkm_tSegType iVolume,
						 char*   isAnnotationFileName,
						 char*   isColorFileName ) { 

  tkm_tErr     eResult                 = tkm_tErr_NoErr;
  Surf_tErr    eSurface                = Surf_tErr_NoErr;
  Volm_tErr    eVolume                 = Volm_tErr_NoErr;
  mriVolumeRef newVolume               = NULL;
  mriVolumeRef newChangedVolume        = NULL;
  MRIS*        mris                    = NULL;
  int          nVertex                 = 0;
  VERTEX*      pVertex                 = 0;
  xColor3n     color;
  int          nStructure              = 0;
  xVoxel       surfRAS;
  xVoxel       anaIdx;
  float        dx                      = 0;
  float        dy                      = 0;
  float        dz                      = 0;
  float        len                     = 0;
  float        d                       = 0;
  char         sError[tkm_knErrStringLen] = "";

  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isAnnotationFileName),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isColorFileName),
		     eResult, tkm_tErr_InvalidParameter );

  /* Make sure we have a surface loaded. */
  DebugAssertThrowX( (NULL != gSurface[tkm_tSurfaceType_Main]),
		     eResult, tkm_tErr_SurfaceNotLoaded );
  
  /* Create the new volume from the existing anatomical. Set it to all
     zeroes.*/
  DebugNote( ("Creating new segmentation volume") );
  eVolume = Volm_New( &newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newVolume, 
				   gAnatomicalVolume[tkm_tVolumeType_Main] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_SetAllValues( newVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Allocate flag volume from the new segmentation volume. Set it to
     all zeroes. */
  DebugNote( ("Creating segmentation flag volume") );
  eVolume = Volm_New( &newChangedVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_CreateFromVolume( newChangedVolume, newVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_CouldntLoadSegmentation );
  eVolume = Volm_SetAllValues( newChangedVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Try to load the color table. */
  DebugNote( ("Loading color table.") );
  eResult = LoadSegmentationColorTable( iVolume, isColorFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eResult),
		     eResult, tkm_tErr_CouldntLoadColorTable );

  /* Read the annotation into the surface. */
  DebugNote( ("Reading annotation.") );
  eSurface = Surf_LoadAnnotation( gSurface[tkm_tSurfaceType_Main],
				  isAnnotationFileName );
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
		     eResult, tkm_tErr_ErrorAccessingSurface );

  /* Get a direct pointer to the surface. */
  DebugNote( ("Getting MRIS") );
  eSurface = Surf_GetMRIS( gSurface[tkm_tSurfaceType_Main], &mris );
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
		     eResult, tkm_tErr_ErrorAccessingSurface );

  /* For every vertex in the surface, get the annot value. Unpack it
     to an RGB value. Look in the LUT to get the corresponding
     index. Get the orig coords from the vertex and set the
     segmentation value at those coords to the index we got. */
  fprintf( stdout, "Converting annotation..." );
  for( nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
    pVertex = &mris->vertices[nVertex];
    
    /* Get the color and then the index. */
    if ( 0 != pVertex->annotation ) {
      MRISAnnotToRGB( pVertex->annotation, 
		      color.mnRed, color.mnGreen, color.mnBlue );
      CLUT_GetIndex( gColorTable[iVolume], &color, &nStructure );
    }

    /* Set all the voxels between the white (current) vertex and the
       pial vertex coords. Find the direction from the pial to the
       normal vertex positions. Then step from the normal outward in
       that direction in 0.1 increments to a distance of 1, convert
       the coord to ana idx, and set the voxel in the seg volume to
       the structure index we found above. */
    dx = pVertex->cx - pVertex->x;
    dy = pVertex->cy - pVertex->y;
    dz = pVertex->cz - pVertex->z;
    
    len = sqrt(dx*dx + dy*dy + dz*dz) ;
    dx /= len; 
    dy /= len; 
    dz /= len;
    
    for( d = 0 ; d <= len; d = d+0.1 ) {
      xVoxl_SetFloat( &surfRAS, 
		      pVertex->x + (d * dx),
		      pVertex->y + (d * dy),
		      pVertex->z + (d * dz) );
      Volm_ConvertRASToIdx( newVolume, &surfRAS, &anaIdx );
      Volm_SetValueAtIdx( newVolume, &anaIdx, (float)nStructure );
    }
    if( !(nVertex % 1000) ) {
      fprintf( stdout, "\rConverting annotation... %.2f%% done",
	       ((float)nVertex / (float)mris->nvertices) * 100.0 );
      fflush( stdout );
    }
  }
  fprintf( stdout, "\rConverting annotation... 100%% done       \n" );


  /* free existing segmentations if present. */
  if( NULL != gSegmentationVolume[iVolume] ) {
    DebugNote( ("Deleting existing roi group") );
    eVolume = Volm_Delete( &gSegmentationVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
		       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
    
    DebugNote( ("Deleting existing flag volume") );
    eVolume = Volm_Delete( &gSegmentationChangedVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
           eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }

  /* Save these volumes. */
  gSegmentationVolume[iVolume]        = newVolume;
  gSegmentationChangedVolume[iVolume] = newChangedVolume;

  /* enable our segmentation options and set segmentation volume in window */
  if( gMeditWindow ) {
    if( tkm_tSegType_Main == iVolume ) {
      DebugNote( ("Enabling segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowSegmentationOptions, "1" );
    } else {
      DebugNote( ("Enabling ayx segmentation options in interface") );
      tkm_SendTclCommand( tkm_tTclCommand_ShowAuxSegmentationOptions, "1" );
    }
  
    DebugNote( ("Setting main segmentation in main window") );
    MWin_SetSegmentationVolume( gMeditWindow, iVolume, 
		      -1, gSegmentationVolume[iVolume] );
    MWin_SetSegmentationColorTable( gMeditWindow, iVolume, 
				    -1, gColorTable[iVolume] );
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr || tkm_tErr_CouldntLoadColorTable,
		   tkm_GetErrorString );
  
  if( eResult != tkm_tErr_CouldntLoadColorTable ) {
  
    xUtil_snprintf( sError, sizeof(sError),
		    "Loading Annotation %s", isAnnotationFileName );
    tkm_DisplayError( sError,
		      tkm_GetErrorString(eResult),
		      "Tkmedit couldn't read the annotation you "
		      "specified. This could be because the annotation "
		      "volume wasn't valid or couldn't be found." );
  }
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;


}

tkm_tErr LoadSegmentationColorTable ( tkm_tSegType iVolume,
				      char*        isColorFileName ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  char      sColorFileName[tkm_knPathLen]   = "";
  CLUT_tErr eColorTable         = CLUT_tErr_NoErr;
  int      nNumEntries         = 0;
  int      nEntry         = 0;
  char      sLabel[CLUT_knLabelLen]     = "";
  char      sTclArguments[tkm_knTclCmdLen]   = "";
  char      sError[tkm_knErrStringLen]     = "";
  
  DebugEnterFunction( ("LoadSegmentationColorTable( iVolume=%d, "
		       "isColorFileName=%s )", iVolume, isColorFileName) );
  
  DebugAssertThrowX( (NULL != isColorFileName),
         eResult, tkm_tErr_InvalidParameter );
  
  /* make a file name */
  DebugNote( ("Making file name from %s", isColorFileName) );
  MakeFileName( isColorFileName, tkm_tFileName_Segmentation, 
    sColorFileName, sizeof(sColorFileName) );
  
  /* try to load color table */
  DebugNote( ("Loading color table") );
  eColorTable = CLUT_New( &gColorTable[iVolume], isColorFileName );
  DebugAssertThrowX( (CLUT_tErr_NoErr == eColorTable),
         eResult, tkm_tErr_CouldntLoadColorTable );
  
  /* build the color table for the interface. go through the entries
     and for each one, make a string of its number and its label, and send
     that as a new entry. */
  DebugNote( ("Clearing color table in interface") );
  tkm_SendTclCommand( tkm_tTclCommand_ClearParcColorTable, "" );
  DebugNote( ("Getting number of color table entries") );
  CLUT_GetNumEntries( gColorTable[iVolume], &nNumEntries );
  for( nEntry = 0; nEntry < nNumEntries; nEntry++ ) {
    DebugNote( ("Getting label for entry %d/%d", nEntry, nNumEntries) );
    eColorTable = CLUT_GetLabel( gColorTable[iVolume], nEntry, sLabel );
    DebugAssertThrowX( (CLUT_tErr_NoErr == eColorTable),
		       eResult, tkm_tErr_Unrecoverable );
    if( strcmp( sLabel, "" ) == 0 )
      xUtil_strncpy( sLabel, "None", sizeof(sLabel) );
    DebugNote( ("Making tcl command") );
    xUtil_snprintf( sTclArguments, sizeof(sTclArguments),
        "%d \"%s\"", nEntry, sLabel );
    tkm_SendTclCommand( tkm_tTclCommand_AddParcColorTableEntry,
      sTclArguments );
  }

  tkm_SendTclCommand( tkm_tTclCommand_UpdateSegmentationColorTable,
		      isColorFileName );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError),
		  "Loading Color Table %s", isColorFileName );
  tkm_DisplayError( sError,
		    tkm_GetErrorString(eResult),
		    "Tkmedit couldn't read the color table you "
		    "specified. This could be because the file wasn't"
		    "valid or wasn't found." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


void SetSegmentationAlpha ( float ifAlpha ) {
  
  char sTclArguments[STRLEN] = "";
  
  if( ifAlpha >= 0 && ifAlpha <= 1.0 ) {

    gfSegmentationAlpha = ifAlpha;

    if( NULL != gMeditWindow ) {
      MWin_SetSegmentationAlpha( gMeditWindow, -1, ifAlpha );
      UpdateAndRedraw ();    
    }

  }
  
  sprintf( sTclArguments, "%f", gfSegmentationAlpha );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateSegmentationVolumeAlpha, 
		      sTclArguments ); 
  
}

void GetSegmentationColorAtVoxel ( tkm_tSegType iVolume, 
				   xVoxelRef    ipVoxel,        
				   xColor3fRef  iBaseColor,
				   xColor3fRef  oColor ) {
  
  int        index        = 0;
  CLUT_tErr  eColorTable  = CLUT_tErr_NoErr;
  xColor3f   roiColor;
  
  /* get the index of this voxel */
  GetSegLabel( iVolume, ipVoxel, &index, NULL );
  
  /* If 0, just use the base color */
  if( 0 == index ) {
    *oColor = *iBaseColor;
    goto cleanup;
  }
  
  /* get the color out of the color map */
  eColorTable = CLUT_GetColorFloat( gColorTable[iVolume], index, &roiColor );
  if( CLUT_tErr_NoErr != eColorTable )
    goto error;
  
  /* blend them */
  oColor->mfRed = (gfSegmentationAlpha * roiColor.mfRed) + 
    (float)((1.0-gfSegmentationAlpha) * iBaseColor->mfRed);
  oColor->mfGreen = (gfSegmentationAlpha * roiColor.mfGreen) + 
    (float)((1.0-gfSegmentationAlpha) * iBaseColor->mfGreen);
  oColor->mfBlue = (gfSegmentationAlpha * roiColor.mfBlue) + 
    (float)((1.0-gfSegmentationAlpha) * iBaseColor->mfBlue);
  
  goto cleanup;
  
 error:
  
  DebugPrint( ( "Error in GetSegmentationColor.\n" ) );
  
  *oColor = *iBaseColor;
  
 cleanup:
  return;
}

void GetSegLabel ( tkm_tSegType iVolume,
		   xVoxelRef    ipVoxel, 
		   int*         onIndex,
		   char*        osLabel ) {
  
  CLUT_tErr   eColorTable  = CLUT_tErr_NoErr;
  Volm_tErr   eVolume      = Volm_tErr_NoErr;
  int         index        = 0;
  float       fValue        = 0;
  
  /* get the voxel at this location */
  eVolume = 
    Volm_GetValueAtIdx( gSegmentationVolume[iVolume], ipVoxel, &fValue );
  if( Volm_tErr_NoErr == eVolume ) {
    index = (int)fValue;
  } else {
    index = 0;
  }
  
  /* get the label out of the table and return it and in the index */
  if( NULL != osLabel ) {
    if( 0 == index ) {
      strcpy( osLabel, "None" );
    } else {
      eColorTable = CLUT_GetLabel( gColorTable[iVolume], index, osLabel );
      if( CLUT_tErr_NoErr != eColorTable ) {
	/* pass an out of bounds notice. */
	strcpy( osLabel, "Out of bounds." );
      }
    }
  }
  
  /* return the index */
  if( NULL != onIndex ) {
    *onIndex = index;
  }
  
  goto cleanup;
  
  goto error;
 error:
  
  DebugPrint( ( "Error in GetSegmentationColor.\n" ) );
  
  if( NULL != onIndex )
    *onIndex = -1;
  if( NULL != osLabel )
    strcpy( osLabel, "None" );
  
 cleanup:
  return;
}

Volm_tVisitCommand AddSimilarVoxelToSelection ( xVoxelRef iAnaIdx,
						float     iValue,
						void*     ipnTarget ) {
  int nIndex = 0;
  int nTargetIndex = 0;
  
  nIndex = (int) iValue;
  nTargetIndex = *(int*)ipnTarget;
  
  if( nIndex == nTargetIndex )
    AddVoxelsToSelection( iAnaIdx, 1 );
  
  return Volm_tVisitComm_Continue;
}

Volm_tVisitCommand AddSimilarVoxelToGraohAvg ( xVoxelRef iAnaIdx,
					       float     iValue,
					       void*     ipnTarget ) {
  int nIndex = 0;
  int nTargetIndex = 0;
  
  nIndex = (int) iValue;
  nTargetIndex = *(int*)ipnTarget;
  
  if( nIndex == nTargetIndex )
    FunV_AddAnatomicalVoxelToSelectionRange( gFunctionalVolume, iAnaIdx );
  
  return Volm_tVisitComm_Continue;
}

Volm_tVisitCommand SetChangedSegmentationValue ( xVoxelRef iAnaIdx,
						 float     iValue,
						 void*     iData ) {
  tkm_tExportSegmentationParamsRef params = NULL;
  float                            value;

  if( iValue ) {

    params = (tkm_tExportSegmentationParamsRef)iData;

    /* Get the value from the segmentation volume and set it in the
       volume to export. */
    Volm_GetValueAtIdx( params->mSrcVolume, iAnaIdx, &value );
    Volm_SetValueAtIdx( params->mDestVolume, iAnaIdx, value );
  }
  
  return Volm_tVisitComm_Continue;
}


tkm_tErr SelectSegLabel ( tkm_tSegType iVolume,
			  int          inIndex ) {
  
  tkm_tErr  eResult     = tkm_tErr_NoErr;
  Volm_tErr eVolume     = Volm_tErr_NoErr;
  int       nNumEntries = 0;
  
  DebugEnterFunction( ("SelectSegLabel ( iVolume=%d, inIndex=%d )", 
		       iVolume, inIndex) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );

  /* make sure we're loaded */
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != gSegmentationVolume || 
		      NULL != gColorTable[iVolume]),
		     eResult, tkm_tErr_SegmentationNotLoaded );
  
  /* check entry index */
  DebugNote( ("Getting number of entries") );
  CLUT_GetNumEntries( gColorTable[iVolume], &nNumEntries );
  DebugAssertThrowX( (inIndex > 0 && inIndex <= nNumEntries),
         eResult, tkm_tErr_InvalidParameter );
  
  /* do it */
  OutputPrint "Selecting... " EndOutputPrint;
  eVolume = 
    Volm_VisitAllVoxels( gSegmentationVolume[iVolume], 
			 &AddSimilarVoxelToSelection, (void*)&inIndex );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  OutputPrint "done.\n" EndOutputPrint;
  
  UpdateAndRedraw ();
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr GraphSegLabel ( tkm_tSegType iVolume,
			 int          inIndex ) {
  
  tkm_tErr  eResult        = tkm_tErr_NoErr;
  Volm_tErr eVolume        = Volm_tErr_NoErr;
  FunV_tErr eFunctional        = FunV_tErr_NoError;
  int      nNumEntries        = 0;
  char      sLabel[CLUT_knLabelLen] = "";
  
  DebugEnterFunction( ("GraphSegLabel ( iVolume=%d, inIndex=%d )", 
		       iVolume, inIndex) );
  
  /* make sure we're loaded */
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != gSegmentationVolume ||
		      NULL != gColorTable[iVolume]),
		     eResult, tkm_tErr_SegmentationNotLoaded );
  
  /* check entry index */
  DebugNote( ("Getting number of entries") );
  CLUT_GetNumEntries( gColorTable[iVolume], &nNumEntries );
  DebugAssertThrowX( (inIndex > 0 && inIndex <= nNumEntries),
         eResult, tkm_tErr_InvalidParameter );
  
  /* begin the selection */
  eFunctional = FunV_BeginSelectionRange( gFunctionalVolume );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
         eResult, tkm_tErr_ErrorAccessingFunctionalVolume );
  
  /* do it */
  OutputPrint "Finding voxels... " EndOutputPrint;
  eVolume = 
    Volm_VisitAllVoxels( gSegmentationVolume[iVolume], 
			 &AddSimilarVoxelToGraohAvg, (void*)&inIndex );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
         eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  OutputPrint "done.\n" EndOutputPrint;
  
  /* finish the selection */
  eFunctional = FunV_EndSelectionRange( gFunctionalVolume );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
         eResult, tkm_tErr_ErrorAccessingFunctionalVolume );
  
  /* set the graph window */
  CLUT_GetLabel( gColorTable[iVolume], inIndex, sLabel );
  if( sLabel != "" )
    FunV_SetLocationString( gFunctionalVolume, sLabel );
  
  UpdateAndRedraw ();
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

static void
RecomputeUpdateCallback( MRI* iMRIValues ) {

  MRIcopy( iMRIValues, 
	   gSegmentationVolume[gRecaluclatingSegemention]->mpMriValues );
  if( gbAcceptingTclCommands && gDisplayIntermediateResults ) {
    UpdateAndRedraw() ;
    MWin_ForceRedraw( gMeditWindow ) ;
  }
}

void RecomputeSegmentation ( tkm_tSegType iVolume ) {
  
  tkm_tErr     eResult   = tkm_tErr_NoErr;
  Volm_tErr    eVolume   = Volm_tErr_NoErr;

  DebugEnterFunction( ("RecomputeSegmentation( iVolume=%d )", 
		       (int)iVolume ) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL == gGCAVolume && NULL == gGCATransform),
		     eResult, tkm_tErr_GCANotLoaded );
  
  /* Delete the old backup if it exists. */
  if( NULL != gPreviousSegmentationVolume[iVolume] ) {
    DebugNote( ("Deleting old backup segmentation") );
    eVolume = Volm_Delete( &gPreviousSegmentationVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), eVolume, 
		       tkm_tErr_ErrorAccessingSegmentationVolume );
  }
  
  /* Make backup volume. */
  DebugNote( ("Creating backup segmentation") );
  eVolume = Volm_DeepClone( gSegmentationVolume[iVolume],
			    &gPreviousSegmentationVolume[iVolume] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume), eVolume, 
		     tkm_tErr_ErrorAccessingSegmentationVolume );
  
  /* Reclassify using the update function declared above. */
  gRecaluclatingSegemention = iVolume;
  GCAreclassifyUsingGibbsPriors( 
             gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues, 
	     gGCAVolume, 
	     gSegmentationVolume[iVolume]->mpMriValues,
	     gGCATransform, 10, 
	     gSegmentationChangedVolume[iVolume]->mpMriValues, 1, 
	     RecomputeUpdateCallback);
  
  UpdateAndRedraw() ;

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  DebugCatchError( eVolume, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

int EditSegmentation ( tkm_tSegType iVolume,
		       xVoxelRef    iAnaIdx,
		       int          inIndex ) {

  int        nOldValue = 0;
  Volm_tValue newValue = 0;
  static int nVolumeIndexBugs = 0;

  DebugEnterFunction( ("EditSegmentation( iVolume=%d, iAnaIdx=%d,%d,%d, "
		       "inIndex=%d )", iVolume, xVoxl_ExpandInt( iAnaIdx ),
		       inIndex) );

  /* For some reason, David was getting random bugs where iVolume
     would be way out of bounds. I'm trying to track this down while
     still letting him work. */
  if( iVolume != tkm_tSegType_Main ) {
    nVolumeIndexBugs++;
    if( nVolumeIndexBugs == 1 ) {
      fprintf( stderr, 
	       "ATTENTION: Please send the .xdebug_tkmedit file\n"
	       "           to Kevin when you're done. Thanks.\n" );
    }
    if( nVolumeIndexBugs < 5 ) {
      xDbg_Printf( "EditSegmentation: iVolume was %d. Stack:\n", iVolume );
      xDbg_PrintStack ();
    }
    iVolume = tkm_tSegType_Main;
  }

  GetSegLabel( iVolume, iAnaIdx, &nOldValue, NULL );
  newValue = (Volm_tValue)inIndex;
  Volm_SetValueAtIdx( gSegmentationVolume[iVolume], iAnaIdx, newValue );
  Volm_SetValueAtIdx( gSegmentationChangedVolume[iVolume], iAnaIdx, 1 );
  
  DebugExitFunction;

  return nOldValue;
}

void CalcSegLabelVolume ( tkm_tSegType iVolume,
			  xVoxelRef    iAnaIdx,
			  int*         onVolume ) {
  
  tkm_tErr  eResult = tkm_tErr_NoErr;
  int      index   = 0;
  tBoolean* visited = NULL;
  int      nSize      = 0;
#ifdef Solaris
  int      i         = 0;
#endif
  
  DebugEnterFunction( ("CalcSegLabelVolume( iAnaIdx=%p, onVolume=%p )",
           iAnaIdx, onVolume) );
  
  DebugAssertThrowX( (NULL != iAnaIdx && NULL != onVolume),
         eResult, tkm_tErr_InvalidParameter );
  
  /* initialize count to 0 */
  *onVolume = 0;
  
  DebugNote( ("Allocating visited tracker volume") );
  nSize = 
    gnAnatomicalDimensionX * 
    gnAnatomicalDimensionY * 
    gnAnatomicalDimensionZ ;
  visited = (tBoolean*) malloc( sizeof(tBoolean) * nSize) ;
  DebugAssertThrowX( (NULL != visited), eResult, tkm_tErr_CouldntAllocate );
  
  /* zero it. solaris doesn't like bzero here... ??? */
#ifdef Solaris
  for( i = 0; i < nSize; i++ )
    visited[i] = FALSE;
#else
  bzero( visited, nSize * sizeof( tBoolean ) );
#endif
  
  DebugNote( ("Getting initial value at %d,%d,%d", 
        xVoxl_ExpandInt( iAnaIdx )) );
  GetSegLabel( iVolume, iAnaIdx, &index, NULL );
  if( 0 == index ) {
    DebugGotoCleanup;
  }
  
  /* recursivly iterate with these values */
  DebugNote( ("Starting iteration") );
  IterateCalcSegLabelVolume( iVolume, iAnaIdx, index, visited, onVolume );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugNote( ("Freeing visited volume") );
  if( NULL != visited )
    free( visited );
  
  DebugExitFunction;
}

void IterateCalcSegLabelVolume ( tkm_tSegType iVolume,
				 xVoxelRef    iAnaIdx, 
				 int          inIndex, 
				 tBoolean*    iVisited, 
				 int*         onVolume ) {
  
  
  Volm_tErr    eVolume   = Volm_tErr_NoErr;
  int         nZ   = 0;
  int         nBeginX   = 0;
  int         nEndX   = 0;
  int         nBeginY   = 0;
  int         nEndY   = 0;
  int         nBeginZ   = 0;
  int         nEndZ   = 0;
  int         nY   = 0;
  int         nX   = 0;
  xVoxel       anaIdx;
  int         index   = 0;
  
  /* if we've already been here, exit */
  if( iVisited[ xVoxl_ExpandToIndex( iAnaIdx,gnAnatomicalDimensionX, 
				     gnAnatomicalDimensionY ) ] == 1 ) {
    goto cleanup;
  }
  
  /* make this voxel as visited */
  iVisited[ xVoxl_ExpandToIndex( iAnaIdx,gnAnatomicalDimensionX, 
				 gnAnatomicalDimensionY ) ] = 1;
  
  
  /* get index at this point */
  GetSegLabel( iVolume, iAnaIdx, &index, NULL );
  
  /* if this is not the voxel we're looking for, exit. */
  if( index != inIndex ) {
    goto cleanup;
  }
  
  /* inc our counter */
  (*onVolume)++;
  
  /* calc our bounds */
  nBeginX = xVoxl_GetX(iAnaIdx) - 1;
  nEndX   = xVoxl_GetX(iAnaIdx) + 1;
  nBeginY = xVoxl_GetY(iAnaIdx) - 1;
  nEndY   = xVoxl_GetY(iAnaIdx) + 1;
  nBeginZ = xVoxl_GetZ(iAnaIdx) - 1;
  nEndZ   = xVoxl_GetZ(iAnaIdx) + 1;
  
  /* check the surrounding voxels */
  for( nZ = nBeginZ; nZ <= nEndZ; nZ++ ) {
    for( nY = nBeginY; nY <= nEndY; nY++ ) {
      for( nX = nBeginX; nX <= nEndX; nX++ ) {
	xVoxl_Set( &anaIdx, nX, nY, nZ );
	eVolume = 
	  Volm_VerifyIdx( gAnatomicalVolume[tkm_tVolumeType_Main], &anaIdx );
	if( Volm_tErr_NoErr == eVolume ) {
	  IterateCalcSegLabelVolume( iVolume, &anaIdx, inIndex, 
				     iVisited, onVolume );
	}
      }
    }
  }
  
 cleanup:
  ;
}

void SetSegmentationValue ( tkm_tSegType iVolume,
			    xVoxelRef    iAnaIdx,
			    int          inIndex ) {
  
  DebugEnterFunction( ("SetSegmentationValue( iVolume=%d, iaAnaIdx=%p "
		       "inIndex=%d )", iVolume, iAnaIdx, inIndex) );

  DebugNote( ("Passing to EditSegmentation") );
  EditSegmentation( iVolume, iAnaIdx, inIndex );

  DebugExitFunction;
}

void SetSegmentationValues ( tkm_tSegType iVolume,
			     xVoxelRef    iaAnaIdx,
			     int          inCount,
			     int          inIndex ) {

  int nVoxel = 0;
  
  DebugEnterFunction( ("SetSegmentationValues( iVolume=%d, iaAnaIdx=%p "
		       "inCount=%d, inIndex=%d )", iVolume, iaAnaIdx,
		       inCount, inIndex) );

  for( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {
    DebugNote( ("EditSegmentation on voxel index %d", inIndex) );
    EditSegmentation( iVolume, &(iaAnaIdx[nVoxel]), inIndex );
  }

  DebugExitFunction;
}


void FloodFillSegmentation ( tkm_tSegType                iVolume,
			     xVoxelRef                   iAnaIdx,
			     tkm_tParcFloodFillParamsRef iParams ) {
  
  tBoolean* visited     = NULL;
  tkm_tErr  eResult     = tkm_tErr_NoErr;
  int       nSize       = 0;
  int       nDimensionX = 0;
  int       nDimensionY = 0;
  int       nDimensionZ = 0;
#ifdef Solaris
  int      i            = 0;
#endif
  
  DebugEnterFunction( ("FloodFillSegmentation( iAnaIdx=%d,%d,%d, iParams=%p )",
		       xVoxl_ExpandInt( iAnaIdx ), iParams) );
  
  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != iAnaIdx),
		     eResult, tkm_tErr_InvalidParameter );

  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume]),
		     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  Volm_GetDimensions( gSegmentationVolume[iVolume], 
		      &nDimensionX, &nDimensionY, &nDimensionZ );
  nSize = nDimensionX * nDimensionY * nDimensionZ;
  
  /* init a volume to keep track of visited voxels */
  DebugNote( ("Allocating visited tracker volume") );
  visited = (tBoolean*) malloc( sizeof(tBoolean) * nSize) ;
  DebugAssertThrowX( (NULL != visited), eResult, tkm_tErr_CouldntAllocate );
  
  /* zero it. solaris doesn't like bzero here... ??? */
#ifdef Solaris
  for( i = 0; i < nSize; i++ )
    visited[i] = FALSE;
#else
  bzero( visited, nSize * sizeof( tBoolean ) );
#endif
  
  /* recursivly iterate with these values */
  DebugNote( ("Starting iteration") );
  IterateFloodFillSegmentation( iVolume, iAnaIdx, iParams, visited );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  UpdateAndRedraw();
  
  DebugNote( ("Freeing visited volume") );
  if( NULL != visited )
    free( visited );
  
  DebugExitFunction;
}

void IterateFloodFillSegmentation ( tkm_tSegType                iVolume,
				    xVoxelRef                   iAnaIdx,
				    tkm_tParcFloodFillParamsRef iParams,
				    tBoolean*                   iVisited ) {
  
  Volm_tErr    eVolume   = Volm_tErr_NoErr;
  int          nZ        = 0;
  int          nBeginX   = 0;
  int          nEndX     = 0;
  int          nBeginY   = 0;
  int          nEndY     = 0;
  int          nBeginZ   = 0;
  int          nEndZ     = 0;
  int          nY        = 0;
  int          nX        = 0;
  xVoxel       anaIdx;
  int          srcValue  = 0;
  float        fValue    = 0;
  float        fDistance = 0;
  
  /* don't do this too many times. */
  iParams->mnIterationDepth++;
  if( iParams->mnIterationDepth > iParams->mnIterationLimit ) {
    iParams->mbIterationLimitReached = TRUE;
    goto cleanup;
  }
  
  /* if we've already been here, exit */
  if( iVisited[ xVoxl_ExpandToIndex( iAnaIdx,gnAnatomicalDimensionX, 
				     gnAnatomicalDimensionY ) ] == 1 ) {
    goto cleanup;
  }
  
  /* make this voxel as visited */
  iVisited[ xVoxl_ExpandToIndex( iAnaIdx,gnAnatomicalDimensionX, 
				 gnAnatomicalDimensionY ) ] = 1;
  
  /* See what src value to get. We can already be sure the proper
     volumes exist because they have been checked by now. */
  switch( iParams->mSource ) {
  case tkm_tVolumeTarget_MainAna:
  case tkm_tVolumeTarget_AuxAna:
    Volm_GetValueAtIdx( gAnatomicalVolume[iParams->mSource], iAnaIdx, 
			&fValue);
    srcValue = (int)fValue;
    break;
  case tkm_tVolumeTarget_MainSeg:
  case tkm_tVolumeTarget_AuxSeg:
    GetSegLabel( iVolume, iAnaIdx, &srcValue, NULL );
    break;
  default:
    return;
    break;
  }

  /* if this is not the voxel we're looking for, exit. */
  if( srcValue < iParams->mnTargetValue - iParams->mnFuzziness ||
      srcValue > iParams->mnTargetValue + iParams->mnFuzziness ) {
    goto cleanup;
  }
  
  /* check distance if >0. if it's over our max distance, exit, */
  if( iParams->mnDistance > 0 ) {
    fDistance = sqrt(
         ((xVoxl_GetX(iAnaIdx) - xVoxl_GetX(&(iParams->mAnchor))) * 
          (xVoxl_GetX(iAnaIdx) - xVoxl_GetX(&(iParams->mAnchor)))) +
         ((xVoxl_GetY(iAnaIdx) - xVoxl_GetY(&(iParams->mAnchor))) * 
          (xVoxl_GetY(iAnaIdx) - xVoxl_GetY(&(iParams->mAnchor)))) +
         ((xVoxl_GetZ(iAnaIdx) - xVoxl_GetZ(&(iParams->mAnchor))) * 
          (xVoxl_GetZ(iAnaIdx) - xVoxl_GetZ(&(iParams->mAnchor)))) );
    if( fDistance > iParams->mnDistance )
      goto cleanup;
  }
  
  /* if so, set this value in the parc */
  SetSegmentationValue( iVolume, iAnaIdx, iParams->mnNewIndex );
  
  /* calc our bounds */
  nBeginX = xVoxl_GetX(iAnaIdx) - 1;
  nEndX    = xVoxl_GetX(iAnaIdx) + 1;
  nBeginY = xVoxl_GetY(iAnaIdx) - 1;
  nEndY    = xVoxl_GetY(iAnaIdx) + 1;
  nBeginZ = xVoxl_GetZ(iAnaIdx) - 1;
  nEndZ    = xVoxl_GetZ(iAnaIdx) + 1;
  
  /* if we're not in 3d, limit one of them based on the orientation */
  if( !iParams->mb3D ) {
    switch( iParams->mOrientation ) {
    case mri_tOrientation_Coronal:
      nBeginZ = nEndZ = xVoxl_GetZ(iAnaIdx);
      break;
    case mri_tOrientation_Horizontal:
      nBeginY = nEndY = xVoxl_GetY(iAnaIdx);
      break;
    case mri_tOrientation_Sagittal:
      nBeginX = nEndX = xVoxl_GetX(iAnaIdx);
      break;
    default:
      goto cleanup;
      break;
    }
  }
  
  /* check the surrounding voxels */
  for( nZ = nBeginZ; nZ <= nEndZ; nZ++ ) {
    for( nY = nBeginY; nY <= nEndY; nY++ ) {
      for( nX = nBeginX; nX <= nEndX; nX++ ) {
	xVoxl_Set( &anaIdx, nX, nY, nZ );
	eVolume = Volm_VerifyIdx( gAnatomicalVolume[tkm_tVolumeType_Main], 
				  &anaIdx );
	if( Volm_tErr_NoErr == eVolume ) {
	  IterateFloodFillSegmentation( iVolume, &anaIdx, iParams, iVisited );
	}
      }
    }
  }
  
 cleanup:
  
  iParams->mnIterationDepth--;
}

void SwapSegmentationAndVolume ( mriVolumeRef   iGroup,
			     mriVolumeRef    iVolume ) {
  
  printf( "UNSUPPORTED\n" );
}

// =============================================================== DTI VOLUMES

tkm_tErr LoadDTIVolume ( char*              isNameEV,
                         char*              isNameFA,
			 tAxis              iRedAxis,
			 tAxis              iGreenAxis,
			 tAxis              iBlueAxis ) {
  
  tkm_tErr     eResult = tkm_tErr_NoErr;
  Volm_tErr    eVolm = Volm_tErr_NoErr;
  char         sError[tkm_knErrStringLen] = "";
  char         sNameWithHint[tkm_knPathLen] = "";
  mriVolumeRef EVVolume = NULL;
  mriVolumeRef FAVolume = NULL;
  xVoxel       screenIdx;
  xVoxel       EVIdx;
  int          zEVX;
  int          zEVY;
  int          zEVZ;
  int          zFAX;
  int          zFAY;
  int          zFAZ;
  float        EVValue;
  float        FAValue;
  int          nFrame;
  
  DebugEnterFunction( ("LoadDTIVolume( isNameEV=%s, isNameFA=%s, iRedAxis=%d, "
		       "iGreenAxis=%d, iBlueAxis=%d )", isNameEV, isNameFA,
		       (int)iRedAxis, (int)iGreenAxis, (int)iBlueAxis ));
  DebugAssertThrowX( (NULL != isNameEV), eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != isNameFA), eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iRedAxis >= 0 && iRedAxis < knNumAxes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iGreenAxis >= 0 && iGreenAxis < knNumAxes),
		     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iBlueAxis >= 0 && iBlueAxis < knNumAxes),
		     eResult, tkm_tErr_InvalidParameter );
  
  /* Create our volume */
  eVolm = Volm_New( &EVVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolm), 
		     eResult, tkm_tErr_ErrorAccessingVolume );
  
  /* First try loading it with just this name. */
  DebugNote( ("Importing data %s", isNameEV) );
  eVolm = Volm_ImportData( EVVolume, isNameEV );
  if( Volm_tErr_NoErr != eVolm ) {
    
    /* Try it again with a BFLOAT hint. */
    DebugNote( ("Copying BFLOAT hint to string") );
    sprintf( sNameWithHint, "%s@BFLOAT", isNameEV );
    
    DebugNote( ("Importing data %s", sNameWithHint) );
    eVolm = Volm_ImportData( EVVolume, sNameWithHint );
    if( Volm_tErr_NoErr != eVolm ) {
      
      /* Try it again with a BSHORT hint. */
      DebugNote( ("Copying BSHORT hint to string") );
      sprintf( sNameWithHint, "%s@BSHORT", isNameEV );
      
      DebugNote( ("Importing data %s", sNameWithHint) );
      eVolm = Volm_ImportData( EVVolume, sNameWithHint );
    }
  }
  
  /* If we still don't have it, throw an error. */
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolm), 
		     eResult, tkm_tErr_CouldntLoadDTIVolume );
  

  /* Now do the same with the FA volume. */
  eVolm = Volm_New( &FAVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolm), 
		     eResult, tkm_tErr_ErrorAccessingVolume );
  
  DebugNote( ("Importing data %s", isNameFA) );
  eVolm = Volm_ImportData( FAVolume, isNameFA );
  if( Volm_tErr_NoErr != eVolm ) {
    
    DebugNote( ("Copying BFLOAT hint to string") );
    sprintf( sNameWithHint, "%s@BFLOAT", isNameFA );
    
    DebugNote( ("Importing data %s", sNameWithHint) );
    eVolm = Volm_ImportData( FAVolume, sNameWithHint );
    if( Volm_tErr_NoErr != eVolm ) {
      
      DebugNote( ("Copying BSHORT hint to string") );
      sprintf( sNameWithHint, "%s@BSHORT", isNameFA );
      
      DebugNote( ("Importing data %s", sNameWithHint) );
      eVolm = Volm_ImportData( FAVolume, sNameWithHint );
    }
  }
  
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolm), 
		     eResult, tkm_tErr_CouldntLoadDTIVolume );

  /* Make sure the dimensions are the same. */
  Volm_GetDimensions( EVVolume, &zEVX, &zEVY, &zEVZ );
  Volm_GetDimensions( FAVolume, &zFAX, &zFAY, &zFAZ );
  DebugAssertThrowX( ((zEVX == zFAX) && (zEVY == zFAY) && (zEVZ == zFAZ)),
		     eResult, tkm_tErr_DTIVolumesDifferentSize );

  /* Now go through the EV volume and assign the processed color value
     scaled by the FA value. (r,g,b) = min(FA,1) * (evx,evy,evz) */
  xVoxl_Set( &EVIdx, 0, 0, 0 );
  do {
    
    /* We're going by the MRI index of the FA volume because it's most
       likely smaller than the screen bounds (256^3) and will be much
       quicker. But we have to convert from MRI idx to screen idx
       before using the index function. Yeah, this is kind of
       backwards. */
    Volm_ConvertMRIIdxToScreenIdx_( EVVolume, &EVIdx, &screenIdx );
      
    Volm_GetValueAtIdxUnsafe( FAVolume, &screenIdx, &FAValue );
    
    for( nFrame = 0; nFrame < 3; nFrame++ ) {
      Volm_GetValueAtIdxFrameUnsafe( EVVolume, &screenIdx, nFrame, &EVValue );
      Volm_SetValueAtIdxFrame( EVVolume, &screenIdx, nFrame, 
			       EVValue * MIN( 1, FAValue) );
    }
    
    if( xVoxl_GetY( &EVIdx ) == 0 ) {
      fprintf( stdout, "\rProcessing DTI volumes: %.2f%% done", 
	       (xVoxl_GetFloatZ( &EVIdx ) / (float)zEVZ) * 100.0 );
    }
    
  } while( xVoxl_IncrementUntilLimits( &EVIdx, zEVX, zEVY, zEVZ ) );
  
  fprintf( stdout, "\rProcessing DTI volumes: 100%% done.         \n" );

  /* Delete the FA volume. */
  Volm_Delete( &FAVolume );

  /* If we already have DTI volumes, delete it. */
  if( NULL != gDTIVolume ) {
    Volm_Delete( &gDTIVolume );
  }
  
  /* Use these DTI volumes. */
  gDTIVolume = EVVolume;
  
  /* Set it in the window. */
  if( NULL != gMeditWindow ) {
    DebugNote( ("Setting DTI volume in main window") );
    MWin_SetDTIVolume( gMeditWindow, -1, gDTIVolume );
  }

  /* Save the axis -> component settings */
  gaDTIAxisForComponent[xColr_tComponent_Red] = iRedAxis;
  gaDTIAxisForComponent[xColr_tComponent_Green] = iGreenAxis;
  gaDTIAxisForComponent[xColr_tComponent_Blue] = iBlueAxis;

  if( NULL != gMeditWindow ) {
    DebugNote( ("Setting DTI axis -> component in main window") );
    MWin_SetDTIAxisForComponent( gMeditWindow, -1, 
				 gaDTIAxisForComponent[xColr_tComponent_Red],
				 xColr_tComponent_Red );
    MWin_SetDTIAxisForComponent( gMeditWindow, -1, 
				 gaDTIAxisForComponent[xColr_tComponent_Green],
				 xColr_tComponent_Green );
    MWin_SetDTIAxisForComponent( gMeditWindow, -1, 
				 gaDTIAxisForComponent[xColr_tComponent_Blue],
				 xColr_tComponent_Blue );
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading DTI volumes %s and %s",
		  isNameEV, isNameFA );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the DTI volume you specified. "
        "This could be because the image format wasn't "
        "recognized, or it couldn't find the proper header, "
        "or the file(s) were unreadable." );
  
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void GetDTIColorAtVoxel ( xVoxelRef        iAnaIdx,
			  mri_tOrientation iPlane,
			  xColor3fRef      iBaseColor,
			  xColor3fRef      oColor ) {
  
  Volm_tErr eVolm = Volm_tErr_NoErr;
  float     x     = 0;
  float     y     = 0;
  float     z     = 0;
  xColr_tComponent comp;
  xColor3f  color;
  
  if( xVoxl_GetX(iAnaIdx) >= 128 &&
      xVoxl_GetZ(iAnaIdx) >= 128 &&
      xVoxl_GetY(iAnaIdx) >= 128 ) {
    ;
  }

  /* Make sure this voxel is in the bounds of the DTI volume */
  eVolm = Volm_VerifyIdxInMRIBounds( gDTIVolume, iAnaIdx );
  if( Volm_tErr_NoErr != eVolm ) {
    *oColor = *iBaseColor;
    return;
  }

  /* Check our 0 alpha special case. */
  if( 0 == gfDTIAlpha ) {
    *oColor = *iBaseColor;
    return;
  }

  /* Get the x, y, and z values. frame 0 has the x value, 1 has y, and
     2 has z. */
  Volm_GetValueAtIdxFrame( gDTIVolume, iAnaIdx, 0, &x );
  Volm_GetValueAtIdxFrame( gDTIVolume, iAnaIdx, 1, &y );
  Volm_GetValueAtIdxFrame( gDTIVolume, iAnaIdx, 2, &z );

  if( x != 0 && y != 0 && z != 0 ) {
    
    /* Map to the colors. We use fabs() here because the value in the
       volume could be negative. */
    for( comp = xColr_tComponent_Red; comp <= xColr_tComponent_Blue; comp++ ) {
      switch( gaDTIAxisForComponent[comp] ) {
      case tAxis_X: xColr_SetFloatComponent( &color, comp, fabs(x) ); break;
      case tAxis_Y: xColr_SetFloatComponent( &color, comp, fabs(y) ); break;
      case tAxis_Z: xColr_SetFloatComponent( &color, comp, fabs(z) ); break;
      default: xColr_SetFloatComponent( &color, comp, 0 ); break;
      }
    }
  
    /* If alpha is 1, just set the color. */
    if( 1 == gfDTIAlpha ) {

      oColor->mfRed = color.mfRed;
      oColor->mfGreen = color.mfGreen;
      oColor->mfBlue = color.mfBlue;

    } else {
      
      /* Blend with the destination color. */
      oColor->mfRed   = MIN( 1, (gfDTIAlpha * color.mfRed) + 
			     (float)((1.0-gfDTIAlpha) * iBaseColor->mfRed));
      oColor->mfGreen = MIN( 1, (gfDTIAlpha * color.mfGreen) + 
			   (float)((1.0-gfDTIAlpha) * iBaseColor->mfGreen));
      oColor->mfBlue  = MIN( 1, (gfDTIAlpha * color.mfBlue) + 
			     (float)((1.0-gfDTIAlpha) * iBaseColor->mfBlue));
    }

  } else {
    *oColor = *iBaseColor;
  }
}
            
void SetDTIAlpha ( float ifAlpha ) {
  
  char sTclArguments[STRLEN] = "";
  
  if( ifAlpha >= 0 && ifAlpha <= 1.0 ) {

    gfDTIAlpha = ifAlpha;

    if( NULL != gMeditWindow ) {
      MWin_SetDTIAlpha( gMeditWindow, -1, ifAlpha );
      UpdateAndRedraw ();    
    }
  }
  
  sprintf( sTclArguments, "%f", gfDTIAlpha );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateDTIVolumeAlpha, 
		      sTclArguments ); 
  
}


// ============================================================== EDITING UNDO

tkm_tErr InitUndoList () {
  
  tkm_tErr   eResult = tkm_tErr_NoErr;
  xUndL_tErr eList   = xUndL_tErr_NoErr;
  
  DebugEnterFunction( ("InitUndoList()") );
  
  /* make the list */
  DebugNote( ("Creating undo list") );
  eList = xUndL_New( &gUndoList, 
         &UndoActionWrapper, &DeleteUndoEntryWrapper );
  DebugAssertThrowX( (xUndL_tErr_NoErr == eList),
         eResult, tkm_tErr_Unrecoverable );
  
  /* set the print function */
  DebugNote( ("Setting print funciton in undo list") );
  eList = xUndL_SetPrintFunction( gUndoList, &PrintEntryWrapper );
  DebugAssertThrowX( (xUndL_tErr_NoErr == eList),
         eResult, tkm_tErr_Unrecoverable );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void DeleteUndoList () {
  
  xUndL_tErr eList = xUndL_tErr_NoErr;
  
  /* delete the list */
  eList = xUndL_Delete( &gUndoList );
  if( xUndL_tErr_NoErr != eList ) {
    DebugPrint( ( "DeleteUndoList(): Error in xUndL_Delete %d: %s\n",
      eList, xUndL_GetErrorString( eList ) ) );
  }
}

void NewUndoEntry ( UndoEntryRef*  opEntry, 
        tUndoActionFunction iFunction,
        xVoxelRef      iVoxel, 
        int        inValue ) {
  
  UndoEntryRef this = NULL;
  
  *opEntry = NULL;
  
  /* allocate the entry */
  this = (UndoEntryRef) malloc( sizeof(UndoEntry) );
  if( NULL == this ) {
    DebugPrint( ( "NewUndoEntry(): Error allocating entry.\n" ) );
    return;
  }
  
  /* vopy the voxel in */
  xVoxl_New( &(this->mVoxel) );
  xVoxl_Copy( this->mVoxel, iVoxel );
  
  /* set the function and value */
  this->mFunction = iFunction;
  this->mnValue    = inValue;
  
  /* return the entry */
  *opEntry = this;
}

void DeleteUndoEntry ( UndoEntryRef* ioEntry ) {
  
  UndoEntryRef this = NULL;
  
  /* see if we have an entry */
  if( NULL == ioEntry ) {
    DebugPrint( ( "DeleteUndoEntry(): NULL parameter.\n" ) );
    return;
  }
  
  this = *ioEntry;
  if( NULL == this ) {
    DebugPrint( ( "DeleteUndoEntry(): Got NULL entry.\n" ) );
    return;
  }
  
  /* delete the voxel */
  xVoxl_Delete( &(this->mVoxel) );
  
  /* delet the entry */
  free( this );
  
  *ioEntry = NULL;
}

void PrintUndoEntry ( UndoEntryRef this ) {
  
  if( NULL == this ) {
    DebugPrint( ( "PrintUndoEntry(): Got NULL entry.\n" ) );
    return;
  }
  
  DebugPrint( ( "%p voxel (%d,%d,%d)  value = %d\n", this,
    xVoxl_ExpandInt(this->mVoxel), this->mnValue ) );
  
}

void ClearUndoList () {
  
  xUndL_tErr eList = xUndL_tErr_NoErr;
  
  /* clear the list */
  eList = xUndL_Clear( gUndoList );
  if( xUndL_tErr_NoErr != eList ) {
    DebugPrint( ( "ClearUndoList(): Error in xUndL_Clear %d: %s\n",
      eList, xUndL_GetErrorString( eList ) ) );
  }
}


void AddVoxelAndValueToUndoList ( tUndoActionFunction iFunction,
          xVoxelRef        iVoxel, 
          int          inValue ) {
  
  UndoEntryRef entry = NULL;
  xUndL_tErr   eList = xUndL_tErr_NoErr;
  
  /* make the entry */
  NewUndoEntry( &entry, iFunction, iVoxel, inValue );
  if( NULL == entry ) {
    DebugPrint( ( "AddVoxelAndValueToUndoList(): Couldn't create entry.\n" ) );
    return;
  }
  
  /* add the entry */
  eList = xUndL_AddEntry( gUndoList, entry );
  if( xUndL_tErr_NoErr != eList ) {
    DebugPrint( ( "AddVoxelAndValueToUndoList(): Error in xUndL_AddEntry "
		  "%d: %s\n", eList, xUndL_GetErrorString( eList ) ) );
  }
}

void RestoreUndoList () {
  
  xUndL_tErr eList = xUndL_tErr_NoErr;
  
  /* restore the list */
  eList = xUndL_Restore( gUndoList );
  if( xUndL_tErr_NoErr != eList ) {
    DebugPrint( ( "RestoreUndoList(): Error in xUndL_Restore %d: %s\n",
      eList, xUndL_GetErrorString( eList ) ) );
  }
  
  /* volume is dirty now */
  SetVolumeDirty( tkm_tVolumeType_Main, TRUE );
  SetVolumeDirty( tkm_tVolumeType_Aux, TRUE );
  
  /* force a redraw in the window */
  UpdateAndRedraw();
}

void DeleteUndoEntryWrapper ( xUndL_tEntryPtr* ipEntryToDelete ) {
  
  if( NULL == *ipEntryToDelete ) {
    DebugPrint( ( "DeleteUndoEntryWrapper(): Got null entry.\n" ) );
    return;
  }
  
  /* typecast and call the delete function */
  DeleteUndoEntry( (UndoEntryRef*)ipEntryToDelete );
  *ipEntryToDelete = NULL;
}

void UndoActionWrapper ( xUndL_tEntryPtr  iUndoneEntry, 
       xUndL_tEntryPtr* opNewEntry ) {
  
  UndoEntryRef entryToUndo = NULL;
  UndoEntryRef undoneEntry = NULL;
  int         value     = 0;
  
  /* we're getting an undo entry that should
     be undone or restored. we also want to
     create a new entry and pass it back, so
     we can undo the undo. */
  
  /* get the entry and check it */
  entryToUndo = (UndoEntryRef)iUndoneEntry;
  if( NULL == entryToUndo ) {
    DebugPrint( ( "UndoActionWrapper(): Got null entry.\n" ) );
    return;
  }
  
  /* call the function for this action. it will return the old value */
  value = entryToUndo->mFunction( entryToUndo->mVoxel, entryToUndo->mnValue );
  
  /* create a new entry, same as the undone one but with the replaced
     value */
  NewUndoEntry( &undoneEntry, 
    entryToUndo->mFunction, entryToUndo->mVoxel, value );
  
  /* pass back the new entry with the replaced value */
  *opNewEntry = undoneEntry;
}

void PrintEntryWrapper ( xUndL_tEntryPtr iEntry ) {
  
  PrintUndoEntry( (UndoEntryRef)iEntry );
}

/* =========================================================== UNDO VOLUME */

UndoVolumeEntryRef NewUndoVolumeEntry ( int iValue ) {
  
  UndoVolumeEntryRef entry = NULL;
  
  entry = (UndoVolumeEntryRef) malloc( sizeof( UndoVolumeEntry ) );
  if( NULL == entry ) {
    DebugPrint( ("NewUndoVolumeEntry: Couldn't create undo volume entry.\n") );
    goto cleanup;
  }
  
  entry->mRestoreValue = iValue;
  
 cleanup:
  
  return entry;
}

void DeleteUndoVolumeEntry ( UndoVolumeEntryRef* iEntry ) {
  
  UndoVolumeEntryRef entry = NULL;
  
  if( NULL == iEntry )
    return;
  
  entry = *iEntry;
  if( NULL == entry )
    return;
  
  free( entry );
  
  *iEntry = NULL;
}

void DeleteUndoVolumeEntryWrapper ( void* iEntry ) {
  
  UndoVolumeEntryRef entry = NULL;
  entry = (UndoVolumeEntryRef)iEntry;
  DeleteUndoVolumeEntry( &entry );
}

tkm_tErr InitUndoVolume () {
  
  tkm_tErr   eResult = tkm_tErr_NoErr;
  xSVol_tErr eVolume = xSVol_tErr_NoErr;
  MRI       *mri_main ;
  
  mri_main = gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ;
  
  DebugEnterFunction( ("InitUndoVolume()") );
  
  /* delete the volume if it exists */
  if( NULL != gUndoVolume )
    DeleteUndoVolume();
  
  /* allocate the volume */
  DebugNote( ("Creating undo volume") );
  eVolume = xSVol_New( &gUndoVolume, mri_main->width, mri_main->height, 
           mri_main->depth );
  DebugAssertThrowX( (xSVol_tErr_NoErr == eVolume),
         eResult, tkm_tErr_Unrecoverable );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
  
}

void DeleteUndoVolume () {
  
  xSVol_tErr eVolume = xSVol_tErr_NoErr;
  
  if( NULL != gUndoVolume ) {
    
    eVolume = xSVol_Delete( &gUndoVolume, DeleteUndoVolumeEntryWrapper );
    if( xSVol_tErr_NoErr != eVolume ) {
      DebugPrint( ( "DeleteUndoVolume: Couldn't delete the undo volume.\n"
        ) );
    }
  }
}

void AddAnaIdxAndValueToUndoVolume ( xVoxelRef    iAnaIdx,
             int    iValue ) {
  
  UndoVolumeEntryRef entry = NULL;
  xSVol_tErr       eVolume = xSVol_tErr_NoErr;
  
  if( NULL == gUndoVolume ) {
    DebugPrint( ("AddAnaIdxAndValueToUndoVolume: Undo volume not inited.\n") );
    goto cleanup;
  }
  
  /* try getting a voxel at this location */
  eVolume = xSVol_Get( gUndoVolume, iAnaIdx, (void**)&entry );
  if( xSVol_tErr_NoErr != eVolume )
    goto cleanup;
  
  /* if it exists, delete it */
  if( NULL != entry ) {
    eVolume = xSVol_Set( gUndoVolume, iAnaIdx, NULL );
    if( xSVol_tErr_NoErr != eVolume )
      goto cleanup;
    
    DeleteUndoVolumeEntry( &entry );
    entry = NULL;
  }
  
  /* create a new voxel */
  entry = NewUndoVolumeEntry( iValue );
  if( NULL == entry )
    goto cleanup;
  
  /* set the voxel at this location */
  eVolume = xSVol_Set( gUndoVolume, iAnaIdx, entry );
  if( xSVol_tErr_NoErr != eVolume )
    goto cleanup;
  
 cleanup:
  return;
}

tBoolean IsAnaIdxInUndoVolume ( xVoxelRef iAnaIdx ) {
  
  tBoolean       bIsInVolume = FALSE;
  xSVol_tErr       eVolume   = xSVol_tErr_NoErr;
  UndoVolumeEntryRef entry   = NULL;
  
  if( NULL == gUndoVolume ) {
    DebugPrint( ( "IsAnaIdxInUndoVolume: Undo volume not inited.\n"
      ) );
    goto cleanup;
  }
  
  /* try getting a voxel at this location */
  eVolume = xSVol_Get( gUndoVolume, iAnaIdx, (void**)&entry );
  if( xSVol_tErr_NoErr != eVolume )
    goto cleanup;
  
  /* if we got on, it's in the volume */
  if( NULL != entry )
    bIsInVolume = TRUE;
  
 cleanup:
  
  return bIsInVolume;
}

void RestoreUndoVolumeAroundAnaIdx ( xVoxelRef iAnaIdx ) {
  
  Volm_tErr       eVolume = Volm_tErr_NoErr;
  UndoVolumeEntryRef entry   = NULL;
  int         nZ       = 0;
  int         nY       = 0;
  int         nX       = 0;
  xVoxel       anaIdx;
  
  /* if this voxel is in the volume... */
  xSVol_Get( gUndoVolume, iAnaIdx, (void**)&entry );
  if( NULL != entry ) {
    
    /* restore the value */
    Volm_SetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Main], iAnaIdx,
      (Volm_tValue)entry->mRestoreValue );
    
    /* remove the voxel */
    xSVol_Set( gUndoVolume, iAnaIdx, NULL );
    DeleteUndoVolumeEntry( &entry );
    
  } else {
    return;
  }
  
  /* try restoring surrounding voxels as well. */
  for( nZ = xVoxl_GetZ(iAnaIdx)-1; nZ <= xVoxl_GetZ(iAnaIdx)+1; nZ++ )
    for( nY = xVoxl_GetY(iAnaIdx)-1; nY <= xVoxl_GetY(iAnaIdx)+1; nY++ )
      for( nX = xVoxl_GetX(iAnaIdx)-1; nX <= xVoxl_GetX(iAnaIdx)+1; nX++ ) {
  xVoxl_Set( &anaIdx, nX, nY, nZ );
  eVolume = Volm_VerifyIdx( gAnatomicalVolume[tkm_tVolumeType_Main], 
          &anaIdx );
  if( Volm_tErr_NoErr == eVolume ) {
    RestoreUndoVolumeAroundAnaIdx( &anaIdx );
  }
      }
  
  UpdateAndRedraw();
}

void ClearUndoVolume () {
  
  if( NULL == gUndoVolume ) {
    DebugPrint( ("ClearUndoVolume: Undo volume not inited.\n") );
    return;
  }
  
  xSVol_Purge( gUndoVolume, DeleteUndoVolumeEntryWrapper );
  
  UpdateAndRedraw();
}

/* ============================================================ HEAD POINTS */

tkm_tErr LoadHeadPts ( char* isHeadPtsFile, 
           char* isTransformFile ) {
  
  tkm_tErr  eResult        = tkm_tErr_NoErr;
  HPtL_tErr eHeadPts        = HPtL_tErr_NoErr;
  char      sHeadPtsFile[tkm_knPathLen]    = "";
  char      sTransformFile[tkm_knPathLen] = "";
  char*      spTransformFileArg      = NULL;
  char      sError[tkm_knErrStringLen]    = "";
  
  DebugEnterFunction( ("LoadHeadPts( isHeadPtsFile=%s, isTransformFile=%s )",
           isHeadPtsFile, isTransformFile) );
  
  /* make filenames */
  DebugNote( ("Making file name from %s", isHeadPtsFile) );
  MakeFileName( isHeadPtsFile, tkm_tFileName_HeadPoints, 
    sHeadPtsFile, sizeof(sHeadPtsFile) );
  if( NULL != isTransformFile ) {
    DebugNote( ("Making file name from %s", isTransformFile) );
    MakeFileName( isTransformFile, tkm_tFileName_HeadPoints, 
      sTransformFile, sizeof(sTransformFile) );
    spTransformFileArg = sTransformFile;
  } else {
    spTransformFileArg = NULL;
  }
  
  DebugNote( ("Creating head points list") );
  eHeadPts = HPtL_New( &gHeadPoints,
           sHeadPtsFile, spTransformFileArg, gIdxToRASTransform );
  DebugAssertThrowX( (HPtL_tErr_NoErr == eHeadPts),
         eResult, tkm_tErr_CouldntLoadHeadPointsList );
  
  DebugNote( ("Setting head point list in window") );
  MWin_SetHeadPointList( gMeditWindow, -1, gHeadPoints );
  
  /* show head points stuff */
  DebugNote( ("Enabling head point options in interface") );
  tkm_SendTclCommand( tkm_tTclCommand_ShowHeadPointLabel, "1" );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading head points %s", 
      isHeadPtsFile );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the head points list you "
        "specified. This could be because the list format "
        "wasn't valid or the file wasn't found, or the "
        "transform file you specified wasn't valid or "
        "found. " );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void RestoreHeadPoints () {
  
  HPtL_tErr eHeadPts  = HPtL_tErr_NoErr;
  
  eHeadPts = HPtL_RestoreTransform( gHeadPoints );
  if ( HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint( ( "HPtL error %d in RestoreHeadPoints: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) ) );
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
    DebugPrint( ( "HPtL error %d in WriteHeadPointsTransform: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) ) );
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
    DebugPrint( ( "HPtL error %d in WriteHeadPointsFile: %s\n",
      eHeadPts, HPtL_GetErrorString( eHeadPts ) ) );
    goto cleanup;
  }
  
  OutputPrint "Wrote head points file.\n" EndOutputPrint;
  
 cleanup:
  return;
}

void RotateHeadPts ( float ifDegrees, 
         char  isAxis ) {
  
  mri_tOrientation orientation;
  
  /* get the orientation */
  MWin_GetOrientation ( gMeditWindow, &orientation );
  
  switch( isAxis ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_Y );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_Y );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_X );
      break;
    default: goto cleanup;
      break;
    }
    break;
  }
  
  UpdateAndRedraw();
  
 cleanup:
  return;
}

void TranslateHeadPts ( float ifDistance, 
      char  isDirection ) {
  
  mri_tOrientation orientation;
  
  /* get the orientation */
  MWin_GetOrientation ( gMeditWindow, &orientation );
  
  switch( isDirection ) {
  case 'x':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Y );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'y':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Z );
      break;
    default: goto cleanup;
      break;
    }
    break;
  case 'z':
    switch( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints,  -ifDistance, tAxis_X );
      break;
    default: goto cleanup;
      break;
    }
    break;
  }
  
  UpdateAndRedraw();
  
 cleanup:
  return;
}

void SetSelectedHeadPointLabel ( char* isNewLabel ) {
  
  HPtL_tHeadPointRef pHeadPoint = NULL;
  MWin_tErr       eWindow  = MWin_tErr_NoErr;
  
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
    DebugPrint( ( "MWin error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, MWin_GetErrorString( eWindow ) ) );
  }
  
 cleanup:
  
  return;
}

void AlignSelectedHeadPointToAnaIdx ( xVoxelRef iAnaIdx ) {
  
  HPtL_tHeadPointRef pHeadPoint = NULL;
  MWin_tErr       eWindow  = MWin_tErr_NoErr;
  HPtL_tErr       eHeadPts  = HPtL_tErr_NoErr;
  
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
    DebugPrint( ( "MWin error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, MWin_GetErrorString( eWindow ) ) );
  }
  
  if(  HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint( ( "HPtL error %d in SetSelectedHeadPointLabel: %s\n",
      eWindow, HPtL_GetErrorString( eHeadPts ) ) );
  }
  
 cleanup:
  
  return;
  
}

/* ===================================================================== GCA */

tkm_tErr LoadGCA ( char* isGCAFileName, char* isTransformFileName ) {
  
  tkm_tErr  eResult        = tkm_tErr_NoErr;
  char      sGCAFileName[tkm_knPathLen]    = "";
  char      sTransformFileName[tkm_knPathLen]  = "";
  GCA*      gca          = NULL;
  TRANSFORM*trans          = NULL;
  char      sError[tkm_knErrStringLen]    = "";
  
  DebugEnterFunction( ("LoadGCA( isGCAFileName=%s, isTransformFileName=%s )",
           isGCAFileName, isTransformFileName) );
  
  /* make filenames */
  DebugNote( ("Making file name from %s", isGCAFileName) );
  MakeFileName( isGCAFileName, tkm_tFileName_GCA, 
    sGCAFileName, sizeof(sGCAFileName) );
  DebugNote( ("Making file name from %s", isTransformFileName) );
  MakeFileName( isTransformFileName, tkm_tFileName_VolumeTransform, 
    sTransformFileName, sizeof(sTransformFileName) );
  
  /* load gca */
  DebugNote( ("Reading GCA with GCAread") );
  gca = GCAread( sGCAFileName );
  DebugAssertThrowX( (NULL != gca), eResult, tkm_tErr_CouldntLoadGCA );
  
  /* load transform */
  DebugNote( ("Loading TRANSFORM with TransformRead") );
  trans = TransformRead( sTransformFileName );
  DebugAssertThrowX( (NULL != trans), eResult, tkm_tErr_CouldntLoadTransform );
  TransformInvert(trans, gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues) ;
  
  /* set globals */
  gGCAVolume  = gca;
  gGCATransform = trans;
  gca    = NULL;
  trans    = NULL;
  
  /* set in window */
  MWin_SetGCA( gMeditWindow, -1, gGCAVolume, gGCATransform );
  tkm_SendTclCommand( tkm_tTclCommand_ShowGCAOptions, "1" );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading GCA %s", 
      isGCAFileName );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the GCA or transform you "
        "specified. This could be because the format "
        "wasn't valid or the file wasn't found." );
  EndDebugCatch;
  
  if( NULL != gca )
    GCAfree( &gca );
#if 0
  GCAhistoScaleImageIntensities(gGCAVolume, 
        gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues) ;
#endif
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr SaveGCA ( char* isGCAFileName ) {
  
  tkm_tErr  eResult        = tkm_tErr_NoErr;
  char      sGCAFileName[tkm_knPathLen]    = "";
  char      sError[tkm_knErrStringLen]    = "";
  
  DebugEnterFunction( ("SaveGCA( isGCAFileName=%s )",
           isGCAFileName) );
  DebugAssertThrowX( (NULL != gGCAVolume), eResult, tkm_tErr_CouldntLoadGCA );
  
  
  /* make filenames */
  DebugNote( ("Making file name from %s", isGCAFileName) );
  MakeFileName( isGCAFileName, tkm_tFileName_GCA, 
    sGCAFileName, sizeof(sGCAFileName) );
  
  /* Save gca */
  DebugNote( ("Reading GCA with GCAread") );
  GCAwrite(gGCAVolume, sGCAFileName );
  
  
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Saving GCA %s", 
      isGCAFileName );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't save the GCA you specified "
        "specified. This could be because the disk was full "
        "or not writeable." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr DeleteGCA () {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("DeleteGCA ()") );
  
  if( NULL != gGCAVolume )
    GCAfree( &gGCAVolume );
  if( NULL != gGCATransform )
    TransformFree( &gGCATransform );
  
  DebugExitFunction;
  
  return eResult;
}

/* ========================================================================= */

/* ===================================================================== VLI */

tkm_tErr LoadVLIs ( char* isFileName1, char* isFileName2 ) {
  
  tkm_tErr  eResult        = tkm_tErr_NoErr;
  char      sFileName1[tkm_knPathLen] = "";
  char      sFileName2[tkm_knPathLen] = "";
  VLI*      vli1 = NULL;
  VLI*      vli2 = NULL;
  char    sError[tkm_knErrStringLen] = "";
  
  DebugEnterFunction( ("ReadVoxelLabels( vli1_name=%s, vli2_name=%s )",
           isFileName1, isFileName2 ) );
  
  /* make filenames */
  DebugNote( ("Making file name from %s", isFileName1) );
  MakeFileName( isFileName1, tkm_tFileName_GCA, 
    sFileName1, sizeof(sFileName1) );
  DebugNote( ("Making file name from %s", isFileName2) );
  MakeFileName( isFileName2, tkm_tFileName_VolumeTransform, 
    sFileName2, sizeof(sFileName2) );
  
  /* load vlis */
  DebugNote( ("Reading VLI1 with VLread") );
  vli1 = VLread( sFileName1 ) ;
  DebugAssertThrowX( (NULL != vli1), eResult, tkm_tErr_CouldntLoadVLI );
  DebugNote( ("Reading VLI2 with VLread") );
  vli2 = VLread( sFileName2 ) ;
  DebugAssertThrowX( (NULL != vli2), eResult, tkm_tErr_CouldntLoadVLI );
  
  /* set globals */
  gVLI1 = vli1 ;
  gVLI2 = vli2 ;
  
  /* set in window */
  MWin_SetVLIs( gMeditWindow, -1, gVLI1, gVLI2, sFileName1, sFileName2 );
  
  /* set tcl menu items */
  tkm_SendTclCommand( tkm_tTclCommand_ShowVLIOptions, (gVLI1!=NULL)?"1":"0" );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  
  xUtil_snprintf( sError, sizeof(sError), "Loading voxel labels %s, %s", 
      sFileName1, sFileName2 );
  tkm_DisplayError( sError,
        tkm_GetErrorString(eResult),
        "Tkmedit couldn't read the VLIs you "
        "specified. This could be because the format "
        "wasn't valid or the file wasn't found." );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr DeleteVLIs () {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("DeleteVLIs ()") );
  
  if( NULL != gVLI1 )
    VLfree( &gVLI1 );
  if( NULL != gVLI2 )
    VLfree( &gVLI2 );
  
  DebugExitFunction;
  
  return eResult;
}


/* ========================================================================= */

/* =================================================================== TIMER */

tkm_tErr StartTimer () {
  
  tkm_tErr    eResult = tkm_tErr_NoErr;
  int      eSystem = 0;
  struct timeval  curTime;
  
  DebugEnterFunction( ("StartTimer()") );
  
  /* mark the current time. */
  DebugNote( ("Clearing timer") );
  timerclear( &curTime );
  DebugNote( ("Getting time with gettimeofday") );
  eSystem = gettimeofday( &curTime, NULL );
  DebugAssertThrowX( (0 == eSystem), eResult, tkm_tErr_GetTimeOfDayFailed );
  gnStartTime = curTime.tv_sec;
  
  /* update status flag */
  gbTimerOn = TRUE;
  tkm_SendTclCommand( tkm_tTclCommand_UpdateTimerStatus, "1" );
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  return eResult;
}

tkm_tErr StopTimer () {
  
  tkm_tErr    eResult = tkm_tErr_NoErr;
  int      eSystem = 0;
  struct timeval  curTime;
  int      nElapsedTime = 0;
  char      sSubjectName[tkm_knNameLen] = "";
  FILE* file = NULL;
  int nCurEntry = 0;
  int nNumEntries = 0;
  tBoolean bFound = FALSE;
  int* anTimes = NULL;
  char** asNames = NULL;
  char* psFileName = NULL;
  
  DebugEnterFunction( ("StartTimer()") );
  
  /* get the current time. */
  DebugNote( ("Clearing timer") );
  timerclear( &curTime );
  DebugNote( ("Getting time with gettimeofday") );
  eSystem = gettimeofday( &curTime, NULL );
  DebugAssertThrowX( (0 == eSystem), eResult, tkm_tErr_GetTimeOfDayFailed );
  
  /* update timer status */
  gbTimerOn = FALSE;
  tkm_SendTclCommand( tkm_tTclCommand_UpdateTimerStatus, "0" );
  
  /* get the elapsed time */
  nElapsedTime = curTime.tv_sec - gnStartTime;
  
  /* get the subject name */
  Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main], 
      sSubjectName, sizeof(sSubjectName) );
  
  /* look for the file name override from the environment variable. */
  psFileName = getenv( ksTkTimerDataFileEnvVar );
  if( NULL != psFileName ) {
    xUtil_strncpy( gsTkTimerFileName, psFileName, sizeof(gsTkTimerFileName) );
  }
  
  /* open file and count how many valid entries in it */
  nNumEntries = 0;
  DebugNote( ("Opening file %s for counting", gsTkTimerFileName) );
  file = fopen( gsTkTimerFileName, "r" );
  if( NULL != file ) {
    while( !feof( file ) ) {
      DebugNote( ("Looking for string and number in file") );
      eSystem = fscanf( file, "%*s %*d" );
      if( EOF != eSystem ) 
  nNumEntries ++;
    }
    DebugNote( ("Closing file") );
    fclose( file );
    file = NULL;
  }
  
  /* allocate storage for the entries */
  DebugNote( ("Allocating storage for %d times", nNumEntries) );
  anTimes = (int*) calloc( nNumEntries, sizeof(int) );
  DebugAssertThrowX( (NULL != anTimes), eResult, tkm_tErr_CouldntAllocate );
  DebugNote( ("Allocating storage for %d names", nNumEntries) );
  asNames = (char**) calloc( nNumEntries, sizeof(char*) );
  DebugAssertThrowX( (NULL != asNames), eResult, tkm_tErr_CouldntAllocate );
  for( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Allocating name string %d", nCurEntry) );
    asNames[nCurEntry] = (char*) malloc( sizeof(char) * tkm_knNameLen );
    DebugAssertThrowX( (NULL != asNames[nCurEntry]),
           eResult, tkm_tErr_CouldntAllocate );
  }
  
  /* read in all the entries from the file */
  nCurEntry = 0;
  DebugNote( ("Opening file %s for reading", gsTkTimerFileName) );
  file = fopen( gsTkTimerFileName, "r" );
  if( NULL != file ) {
    while( !feof( file ) ) {
      DebugNote( ("Getting line %d", nCurEntry) );
      eSystem = fscanf( file, "%s %d", 
      asNames[nCurEntry], &(anTimes[nCurEntry]) );
      if( 2 == eSystem ) {
  nCurEntry ++;
      }
    }
    DebugNote( ("Closing file") );
    fclose( file );
    file = NULL;
  }
  
  /* look for the entry with this subject name. if we find it, add the 
     elapsed time to the time entry. */
  bFound = FALSE;
  for( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Looking for %s in entry %d", sSubjectName, nCurEntry) );
    if( strcmp( asNames[nCurEntry], sSubjectName ) == 0 ) {
      anTimes[nCurEntry] += nElapsedTime;
      bFound = TRUE;
      break;
    }
  }
  
  /* now open the file and write out all the entries. if the entry
     was not found, add a new entry at the end with the elapsed time. */
  DebugNote( ("Opening file %s for writing", gsTkTimerFileName) );
  file = fopen( gsTkTimerFileName, "w" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_CouldntWriteFile );
  for( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Writing entry %d to file", nCurEntry) );
    fprintf( file, "%s %d\n", asNames[nCurEntry], anTimes[nCurEntry] );
  }
  if( !bFound ) {
    DebugNote( ("Writing new entry to file") );
    fprintf( file, "%s %d\n", sSubjectName, nElapsedTime );
  }
  DebugNote( ("Closing file") );
  fclose( file );
  file = NULL;
  
  
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  tkm_DisplayError( "Writing TkTimer data",
        "Couldn't write file",
        "Tkmedit couldn't write to the TkTimer data file. "
        "If you specified a file with the TKTIMER_FILENAME "
        "environment variable, please make sure it is valid "
        "and writeable. If you didn't, please make sure that "
        "the directory you loaded Tkmedit from is writeable, "
        "since that's where the file will be made." );
  EndDebugCatch;
  
  if( NULL != asNames ) {
    for( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
      DebugNote( ("Freeing name %d", nCurEntry) );
      free( asNames[nCurEntry] );
    }
    DebugNote( ("Freeing name array") );
    free( asNames );
  }
  if( NULL != anTimes ) {
    DebugNote( ("Freeing times array") );
    free( anTimes );
  }
  if( NULL != file ) {
    DebugNote( ("Closing file") );
    fclose( file );
  }
  
  DebugExitFunction;
  return eResult;
}

/* ============================================================ ACCESS FUNCS */

void tkm_MakeProgressBar ( char* isName, char* isDesc ) {
  
  fprintf( stdout, "%s: %s\n", isName, isDesc );
  fflush( stdout );
}

void tkm_UpdateProgressBar ( char* isName, float ifPercent ) {
  
  fprintf( stdout, "\r%s: %.2f%% done        ", isName, ifPercent );
  fflush( stdout );
}

void tkm_FinishProgressBar ( char* isName ) {
  
  fprintf( stdout, "\r%s: 100%% done       \n", isName );
  fflush( stdout );
}

void tkm_DisplayMessage ( char* isMessage ) {
  
  fprintf( stdout, isMessage );
  fflush( stdout );
}

void tkm_DisplayError ( char* isAction, char* isError, char* isDesc ) {
  
  char sTclArguments[tkm_knTclCmdLen];
  
  DebugEnterFunction( ("tkm_DisplayError( isAction=%s, isError=%s, "
           "isDesc=%s )", isAction, isError, isDesc) );
  
  DebugNote( ("Making error string") );
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), 
      "\"%s\" \"%s\" \"%s\"", isAction, isError, isDesc );
  
  DebugPrint( ("ERROR: %s\n   %s\n", isAction, isError) );
  tkm_SendTclCommand( tkm_tTclCommand_FormattedErrorDlog, sTclArguments );
  
  DebugExitFunction;
}

void tkm_DisplayAlert ( char* isAction, char* isMsg, char* isDesc ) {
  
  fprintf( stdout, "ALERT: %s\n        %s\n   %s\n", 
     isAction, isMsg, isDesc );
  fflush( stdout );
  
  DebugPrint( ("ALERT: %s\n%s\n", isAction, isMsg) );
}

void tkm_MakeControlPoint ( xVoxelRef iAnaIdx ) {
  
  if( NULL == iAnaIdx ) {
    DebugPrint( ( "tkm_NewCtrlPt(): Passed NULL voxel.\n" ) );
    return;
  }
  
  NewControlPoint ( iAnaIdx, TRUE );
}

void tkm_RemoveControlPointWithinDist ( xVoxelRef   iAnaIdx,
          mri_tOrientation iPlane,
          int       inDistance ) {
  
  float         fDistance = 0;
  xVoxelRef    pCtrlPt   = NULL;
  
  /* find the closest control point */
  fDistance = FindNearestControlPoint( iAnaIdx, iPlane, &pCtrlPt );
  
  /* if we found one and it's in range... */
  if( NULL != pCtrlPt &&
      fDistance >= 0.0 &&
      fDistance <= (float)inDistance ) {
    
    /* delete it */
    DeleteControlPoint( pCtrlPt );
  }
}

void tkm_WriteControlFile () {
  
  WriteControlPointFile( );
}

void tkm_EditAnatomicalVolumeInRange( tkm_tVolumeType  iVolume, 
				      xVoxelRef        inVolumeVox, 
				      tVolumeValue     inLow, 
				      tVolumeValue     inHigh, 
				      tVolumeValue     inNewValue ) {
  
  EditAnatomicalVolumeInRangeArray( iVolume, inVolumeVox, 1,
				    inLow, inHigh, inNewValue );
  
}

void tkm_EditAnatomicalVolumeInRangeArray( tkm_tVolumeType  iVolume, 
					   xVoxelRef        iaVolumeVox, 
					   int              inCount,
					   tVolumeValue     inLow, 
					   tVolumeValue     inHigh, 
					   tVolumeValue     inNewValue ) {
  
  EditAnatomicalVolumeInRangeArray( iVolume, iaVolumeVox, inCount,
				    inLow, inHigh, inNewValue );
  
}

void tkm_SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
				     int             iAnaX0,
				     int             iAnaX1,
				     int             iAnaY0,
				     int             iAnaY1,
				     int             iAnaZ0,
				     int             iAnaZ1,
				     float           iNewValue ) {

  SetAnatomicalVolumeRegion( iVolume, 
			     iAnaX0, iAnaX1, 
			     iAnaY0, iAnaY1, 
			     iAnaZ0, iAnaZ1, 
			     iNewValue );
}

void tkm_SelectVoxel ( xVoxelRef iAnaIdx ) {
  
  AddVoxelsToSelection ( iAnaIdx, 1 );
}

void tkm_SelectVoxelArray ( xVoxelRef iaAnaIdx, int inCount ) {
  
  AddVoxelsToSelection ( iaAnaIdx, inCount );
}

void tkm_DeselectVoxel ( xVoxelRef iAnaIdx ) {
  
  RemoveVoxelsFromSelection( iAnaIdx, 1 );
}

void tkm_DeselectVoxelArray ( xVoxelRef iaAnaIdx, int inCount ) {
  
  RemoveVoxelsFromSelection( iaAnaIdx, inCount );
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

void tkm_ClearUndoVolume () {
  
  ClearUndoVolume();
}

void tkm_RestoreUndoVolumeAroundAnaIdx ( xVoxelRef iAnaIdx ) {
  
  RestoreUndoVolumeAroundAnaIdx( iAnaIdx );
}

tBoolean tkm_IsAnaIdxInUndoVolume ( xVoxelRef iAnaIdx ) {
  
  return IsAnaIdxInUndoVolume( iAnaIdx );
}

void tkm_GetHeadPoint ( xVoxelRef      iAnaIdx,
      mri_tOrientation    iOrientation,
      tBoolean      ibFlat,
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
    
    if( ibFlat ) {
      eHeadPts = HPtL_FindFlattenedNearestPoint( gHeadPoints, plane,
             iAnaIdx, &pPoint );
    } else {
      eHeadPts = HPtL_FindNearestPoint( gHeadPoints, plane,
          1.0, iAnaIdx, &pPoint );
    }
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
  
  GotoSavedCursor();
}

void tkm_GetAnaDimension  ( tkm_tVolumeType iVolume,
          int*      onDimensionX, int* onDimensionY, int *onDimensionZ  ) {
  
  Volm_GetDimensions( gAnatomicalVolume[iVolume], onDimensionX, onDimensionY, onDimensionZ );
}

tBoolean tkm_IsValidAnaIdx ( tkm_tVolumeType iVolume,
			     xVoxelRef     iAnaIdx ) {
  
  Volm_tErr    eVolume   = Volm_tErr_NoErr;
  eVolume = Volm_VerifyIdx( gAnatomicalVolume[iVolume], iAnaIdx );
  return ( Volm_tErr_NoErr == eVolume );
}

void tkm_GetSegmentationColorAtVoxel ( tkm_tSegType iVolume,
			      xVoxelRef    inVoxel, 
			      xColor3fRef  iBaseColor,
			      xColor3fRef  oColor ) {
  
  GetSegmentationColorAtVoxel( iVolume, inVoxel, iBaseColor, oColor );
}

void tkm_GetSegLabel ( tkm_tSegType iVolume,
		       xVoxelRef    inVoxel, 
		       int*         onIndex,
		       char*        osLabel ) {
  
  GetSegLabel( iVolume, inVoxel, onIndex, osLabel );
}

void tkm_SelectCurrentSegLabel ( tkm_tSegType iVolume,
				 int          inIndex ) {

  tkm_tErr eResult = tkm_tErr_NoErr;
  
  eResult = SelectSegLabel( iVolume, inIndex );
  if( tkm_tErr_NoErr != eResult ) {
    tkm_DisplayError( "Selecting Current Segmentation Label",
		      "Tool failed",
		      "Tkmedit couldn't select the current label. You are "
		      "probably trying to select in invalid label. Make sure "
		      "you click on a label first and its name appears in "
		      "the tools window in the seg label line." );
  }
}

void tkm_GraphCurrentSegLabelAvg ( tkm_tSegType iVolume,
			      int          inIndex ) {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  eResult = GraphSegLabel( iVolume, inIndex );
  if( tkm_tErr_NoErr != eResult ) {
    tkm_DisplayError( "Graphing Current Segmentation Label Average",
		      "Tool failed",
		      "Tkmedit couldn't graph the current label. You are "
		      "probably trying to graph in invalid label. Make sure "
		      "you click on a label first and its name appears in "
		      "the tools window in the seg label line. Also, make "
		      "sure that functional time course data is loaded." );
  }
}

void tkm_CalcSegLabelVolume ( tkm_tSegType iVolume,
			      xVoxelRef    iAnaIdx,
			      int*         onVolume ) {
  
  CalcSegLabelVolume( iVolume, iAnaIdx, onVolume );
}

void tkm_EditSegmentation ( tkm_tSegType iVolume,
			    xVoxelRef    iAnaIdx,
			    int          inIndex ) {
  
  DebugEnterFunction( ("tkm_EditSegmentation( iVolume=%d, iAnaIdx=%p, "
		       "inIndex=%d )", iVolume, iAnaIdx, inIndex) );

  SetSegmentationValue( iVolume, iAnaIdx, inIndex );

  DebugExitFunction;
}

void tkm_EditSegmentationArray ( tkm_tSegType iVolume,
				 xVoxelRef    iaAnaIdx,
				 int          inCount,
				 int          inIndex ) {
  
  DebugEnterFunction( ("tkm_EditSegmentationArray( iVolume=%d, iaAnaIdx=%p "
		       "inCount=%d, inIndex=%d )", iVolume, iaAnaIdx,
		       inCount, inIndex) );

  SetSegmentationValues( iVolume, iaAnaIdx, inCount, inIndex );

  DebugExitFunction;
}


void tkm_FloodFillSegmentation ( tkm_tSegType    iVolume,
				 xVoxelRef       iAnaIdx,
				 int             inIndex,
				 tBoolean        ib3D,
				 tkm_tVolumeType iSrc,
				 int             inFuzzy,
				 int             inDistance ) {
  
  tkm_tParcFloodFillParams params;
  char                     sTclArguments[1024] = "";
  float                    anaValue;
  
  params.mnFuzziness             = inFuzzy;
  params.mnDistance              = inDistance;
  xVoxl_Copy( &params.mAnchor, iAnaIdx );
  params.mb3D                    = ib3D;
  params.mSource                 = iSrc;
  params.mnNewIndex              = inIndex;
  params.mnIterationDepth        = 0;
  params.mnIterationLimit        = knMaxFloodFillIteration;
  params.mbIterationLimitReached = FALSE;
  MWin_GetOrientation ( gMeditWindow, &params.mOrientation );
  
  /* see what target value to get. for everything other than the main
     anatomical volume, check to see if it exists first. for the
     segmentation targets, set the fuzzinees to 0, since it doesn't
     really apply to label values. */
  switch( params.mSource ) {
  case tkm_tVolumeTarget_MainAna:
    Volm_GetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Main], 
			iAnaIdx, &anaValue );
    params.mnTargetValue = (int)anaValue;
    break;
  case tkm_tVolumeTarget_AuxAna:
    if( NULL == gAnatomicalVolume[tkm_tVolumeType_Aux] ) {
      strcpy( sTclArguments, 
	      "\"Cannot use aux volume as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    Volm_GetValueAtIdx( gAnatomicalVolume[tkm_tVolumeType_Aux], 
			iAnaIdx, &anaValue);
    params.mnTargetValue = (int)anaValue;
    break;
  case tkm_tVolumeTarget_MainSeg:
    if( NULL == gSegmentationVolume[tkm_tSegType_Main] ) {
      strcpy( sTclArguments, 
	      "\"Cannot use segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    GetSegLabel( tkm_tSegType_Main, iAnaIdx, &params.mnTargetValue, NULL );
    params.mnFuzziness = 0;
    break;
  case tkm_tVolumeTarget_AuxSeg:
    if( NULL == gSegmentationVolume[tkm_tSegType_Aux] ) {
      strcpy( sTclArguments, 
	    "\"Cannot use aux segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    GetSegLabel( tkm_tSegType_Aux, iAnaIdx, &params.mnTargetValue, NULL );
    params.mnFuzziness = 0;
    break;
  default:
    return;
    break;
  }
  
  /* do it! */
  FloodFillSegmentation( iVolume, iAnaIdx, &params );
  
  /* see if our iteration limit was reached. if so, spit out an error. */
  if( params.mbIterationLimitReached ) {
    DebugPrint( ( "FloodFillSegmentation: Recursion limit reached.\n" ) );
    strcpy( sTclArguments, 
      "\"Flood fill recursion limit reached. The region you wished to "
      "fill was too big for me to handle. I have filled as much as I "
      "can. You should try filling again from another starting point, "
      "close to the edge of the part that was filled.\nIf you find this "
      "happening a lot, please email the programmer and bug him to "
      "change this.\"" );
    tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
    goto cleanup;
  }
  
 cleanup:
  return;
}

void tkm_GetDTIColorAtVoxel ( xVoxelRef        iAnaIdx,
			      mri_tOrientation iPlane,
			      xColor3fRef      iBaseColor,
            xColor3fRef      oColor ) {
  
  GetDTIColorAtVoxel( iAnaIdx, iPlane, iBaseColor, oColor );
}
            
void tkm_SetSurfaceDistance    ( xVoxelRef iAnaIdx,
         float     ifDistance ) {
  
  if( NULL == gSurface[tkm_tSurfaceType_Main] ) {
    return;
  }
  
  Surf_SetVertexValue( gSurface[tkm_tSurfaceType_Main], 
           Surf_tVertexSet_Main, Surf_tValueSet_Val,
           iAnaIdx, ifDistance );
}

void tkm_ShowNearestSurfaceVertex ( Surf_tVertexSet iVertexSet ) {
  
  if( iVertexSet < 0 || iVertexSet >= Surf_knNumVertexSets ) {
    return;
  }

  FindNearestSurfaceVertex( iVertexSet );
}


char *kTclCommands [tkm_knNumTclCommands] = {
  
  "UpdateLinkedCursorFlag",
  "UpdateVolumeCursor",
  "UpdateVolumeSlice",
  "UpdateRASCursor",
  "UpdateTalCursor",
  "UpdateScannerCursor",
  "UpdateMNICursor",
  "UpdateVolumeName",
  "UpdateVolumeValue",
  "UpdateAuxVolumeName",
  "UpdateAuxVolumeValue",
  "UpdateFunctionalCoords",
  "UpdateFunctionalRASCoords",
  "UpdateFunctionalValue",
  "UpdateSegLabel",
  "UpdateAuxSegLabel",
  "UpdateHeadPointLabel",
  "UpdateSurfaceDistance",
  "UpdateZoomLevel",
  "UpdateOrientation",
  "UpdateDisplayFlag",
  "UpdateTool",
  "UpdateBrushTarget",
  "UpdateBrushShape",
  "UpdateBrushInfo",
  "UpdateCursorColor",
  "UpdateCursorShape",
  "UpdateSurfaceLineWidth",
  "UpdateSurfaceLineColor",
  "UpdateParcBrushInfo",
  "UpdateVolumeColorScaleInfo",
  "UpdateSegmentationVolumeAlpha",
  "UpdateDTIVolumeAlpha",
  "UpdateTimerStatus",
  "UpdateSubjectDirectory",
  "UpdateSegmentationColorTable",
  "UpdateVolumeDirty",
  "UpdateAuxVolumeDirty",
  
  /* display status */
  "ShowVolumeCoords",
  "ShowRASCoords",
  "ShowTalCoords",
  "ShowAuxValue",
  "ShowSegLabel",
  "ShowAuxSegLabel",
  "ShowHeadPointLabel",
  "ShowFuncCoords",
  "ShowFuncValue",
  "tkm_SetMenuItemGroupStatus tMenuGroup_AuxVolumeOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_DirtyAnatomicalVolume",
  "tkm_SetMenuItemGroupStatus tMenuGroup_DirtyAuxAnatomicalVolume",
  "tkm_SetMenuItemGroupStatus tMenuGroup_VolumeMainTransformLoadedOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_VolumeAuxTransformLoadedOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_OverlayOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_TimeCourseOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_SurfaceLoading",
  "tkm_SetMenuItemGroupStatus tMenuGroup_SurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_OriginalSurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_PialSurfaceViewing",
  "tkm_SetMenuItemGroupStatus tMenuGroup_HeadPoints",
  "tkm_SetMenuItemGroupStatus tMenuGroup_VLIOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_GCAOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_DTIOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_Registration",
  "tkm_SetMenuItemGroupStatus tMenuGroup_SegmentationOptions",
  "tkm_SetMenuItemGroupStatus tMenuGroup_AuxSegmentationOptions",
  "ClearParcColorTable",
  "AddParcColorTableEntry",
  
  /* histogram */
  "BarChart_Draw",
  
  /* interface configuration */
  "wm geometry . ",
  "wm deiconify .; raise .",
  "CsurfInterface",
  "tkm_Finish",
  
  /* misc */
  "ErrorDlog",
  "FormattedErrorDlog",
  "AlertDlog",
  "tkm_MakeProgressDlog",
  "tkm_UpdateProgressDlog",
  "tkm_DestroyProgressDlog"
};

void tkm_SendTclCommand ( tkm_tTclCommand inCommand,
        char* inArguments ) {
  
  char theCommand[1024];
  
  if ( inCommand < 0
       || inCommand >= tkm_knNumTclCommands )
    return;
  
  sprintf ( theCommand, "%s %s", kTclCommands[inCommand], inArguments );
  //  DebugPrint( ( "[] %s\n", theCommand ) );
  SendTCLCommand ( theCommand );
}

char* tkm_GetErrorString ( tkm_tErr ieCode ) {
  
  tkm_tErr eCode = ieCode;
  
  if ( ieCode <= tkm_tErr_InvalidErrorCode ||
       ieCode >= tkm_knNumErrorCodes ) {
    eCode = tkm_tErr_InvalidErrorCode;
  }
  
  return tkm_ksaErrorStrings [eCode];
}

void PrintCachedTclErrorDlogsToShell () {
  
#ifndef Solaris
#ifndef IRIX
  int nCommand = 0;
  char sSpace[3] = " ";
  char sCommand[tkm_knTclCmdLen] = "";
  char* pCommand = NULL;
  
  
  /* send all our cached commands */
  for( nCommand = 0; nCommand < gNumCachedCommands; nCommand++ ) {
    
    /* if this is an error command */
    if( strstr( gCachedTclCommands[nCommand], 
    kTclCommands[tkm_tTclCommand_ErrorDlog] ) != NULL ) {
      
      /* get the command, go to the first space, copy that into a string,
   and print it to the shell. */
      pCommand = gCachedTclCommands[nCommand];
      strsep( &pCommand, sSpace );
      xUtil_snprintf( sCommand, sizeof(sCommand), "%s\n", pCommand );
      fprintf( stdout, sCommand );
      fflush( stdout );
    }
  }
#endif
#endif
}

tkm_tErr EnqueueScript ( char* isFileName ) {
  
  tkm_tErr eResult = tkm_tErr_NoErr;
  
  DebugEnterFunction( ("EnqueueScript ( isFileName=%s )", isFileName) );
  
  DebugNote( ("Checking paramaters") );
  DebugAssertThrowX( (isFileName != NULL), 
         eResult, tkm_tErr_InvalidParameter );
  
  DebugAssertThrowX( (gnNumCachedScripts < tkm_knMaxNumCachedScripts),
         eResult, tkm_tErr_CouldntCacheScriptName );
  
  /* cache the script */
  xUtil_strncpy( gsCachedScriptNames[gnNumCachedScripts], isFileName,
     sizeof( gsCachedScriptNames[gnNumCachedScripts] ) );
  gnNumCachedScripts++;
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

tkm_tErr ExecuteQueuedScripts () {
  
  tkm_tErr eResult          = tkm_tErr_NoErr;
  int     eTcl            = TCL_OK;
  int     nScript          = 0;
  char     sError[tkm_knErrStringLen] = "";
  
  DebugEnterFunction( ("ExecuteQueuedScripts ()") );
  
  for( nScript = 0; nScript < gnNumCachedScripts; nScript++ ) {
    
    DebugAssertThrowX( (gsCachedScriptNames[nScript] != NULL),
           eResult, tkm_tErr_InvalidScriptName );
    
    DebugNote( ("Parsing script file %s", gsCachedScriptNames[nScript]) );
    eTcl = Tcl_EvalFile( interp, gsCachedScriptNames[nScript] );
    if( eTcl != TCL_OK ) {
      xUtil_snprintf( sError, sizeof(sError), "%s", interp->result );
      tkm_DisplayError( "Parsing script file",
      sError,
      "This error was encountered while running the "
      "script specified with the -tcl option." );
    } else {
      OutputPrint "Executed script file %s successfully.\n",
  gsCachedScriptNames[nScript] EndOutputPrint;
    }
  }
  
  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

void HandleSegfault ( int nSignal ) {
  
  tBoolean bDirty = FALSE;
  
  printf( "\n\n===================================\n" );
  printf( "ERROR: A segfault has occurred. This is not your fault,\n" );
  printf( "  : but is most likely an unrecoverable error and has\n" );
  printf( "  : made the program unstable.\n" );
  printf( "  :\n" );
  
  IsVolumeDirty( tkm_tVolumeType_Main, &bDirty );
  if( bDirty ) {
    printf( "    : Attempting to save main volume to /tmp directory...\n" );
    SaveVolume( tkm_tVolumeType_Main, "/tmp" );
    printf( "    : Data was saved to /tmp.\n" );
    printf( "    :\n" );
  }
  
  printf( "  : Please send the contents of the file .xdebug_tkmedit\n" );
  printf( "  : that should be in this directory to analysis-bugs@nmr.mgh.harvard.edu\n");
  printf( "  :\n" );
  
  printf( "  : Now exiting...\n" );
  printf( "  :\n" );
  exit( 1 );
}



xList_tCompare CompareVoxels ( void* inVoxelA, void* inVoxelB ) {
  
  xList_tCompare eResult = xList_tCompare_GreaterThan;
  xVoxelRef   voxelA   = NULL;
  xVoxelRef   voxelB   = NULL;
  
  voxelA = (xVoxelRef) inVoxelA;
  voxelB = (xVoxelRef) inVoxelB;
  
  if( xVoxl_IsEqualInt( voxelA, voxelB ) ) {
    eResult = xList_tCompare_Match;
  } else {
    eResult = xList_tCompare_GreaterThan;
  }
  
  return eResult;
}

tkm_tErr
LoadGCARenormalization(char *renormalization_fname)
{
  FILE   *fp ;
  int   *labels, nlines, i ;
  float   *intensities, f1, f2 ;
  char   *cp, line[STRLEN] ;
  
  if (!gGCAVolume)
    return(tkm_tErr_NoErr) ;
  
  
  fp = fopen(renormalization_fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s",
        Progname, renormalization_fname) ;
  
  cp = fgetl(line, 199, fp) ;
  nlines = 0 ;
  while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
  rewind(fp) ;
  printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
  labels = (int *)calloc(nlines, sizeof(int)) ;
  intensities = (float *)calloc(nlines, sizeof(float)) ;
  cp = fgetl(line, 199, fp) ;
  for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ; intensities[i] = f2 ;
      cp = fgetl(line, 199, fp) ;
    }
  GCArenormalizeIntensities(gGCAVolume, labels, intensities, nlines) ;
  free(labels) ; free(intensities) ;
  return(tkm_tErr_NoErr) ;
}
