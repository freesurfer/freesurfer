/**
 * @file  tkmedit.c
 * @brief Tcl/Tk-based MRI volume and surface viewer and editor
 *
 * TkMedit displays anatomical data and allows the user to navigate through
 * that data and view it from different orientations. TkMedit also displays
 * other data types such as functional data and surfaces as overlays onto
 * this anatomical data.
 * See: http://surfer.nmr.mgh.harvard.edu/fswiki/TkMeditGuide
 */
/*
 * Original Author: Martin Sereno and Anders Dale, 1996
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/10/15 15:43:39 $
 *    $Revision: 1.346 $
 *
 * Copyright (C) 2002-2011, CorTechs Labs, Inc. (La Jolla, CA) and
 * The General Hospital Corporation (Boston, MA).
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
 * in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
 * and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#undef VERSION

char *VERSION = "$Revision: 1.346 $";

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
#include <libgen.h>
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "const.h"
#include "mrisurf.h"
#include "transform.h"
#include "version.h"
#include "ctrpoints.h"
#include "tiffio.h"
#include "fio.h"
#include "tkmedit.h"
#include "fsgdf_wrap.h"
#include "fsgdf.h"
#include "mri2.h"

#include <tcl.h>
//#include <tclDecls.h>
#include <tk.h>
//
// It seems that the later version of Tix uses ITcl and ITk
// for RedHat Enterprise Linux only.   Do the following only
// for RedHat Enterprise Linux.
// stupid tix people who cannot handle version.   I had to use gcc version.
// you cannot include itk.h either(producing so many unknowns) either.
#if NEEDS_ITCL_ITK
#ifndef Itcl_Init
int Itcl_Init(Tcl_Interp* interp);
#endif
#ifndef Itk_Init
int Itk_Init(Tcl_Interp* interp);
#endif
#endif

//////////////////////////////////////////////////////
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
  "The environment variable FREESURFER_HOME is not defined.",
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
  "Couldn't load the surface annotationx[.",
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
  "Couldn't load the FSGDF file.",
  "Error accessing a file.",
  "Error accessing the anatomical volume.",
  "Error accessing the RAS transform.",
  "Error accessing the segmentation.",
  "Error accessing the functional volume.",
  "Error accessing the list.",
  "Error accessing surface.",
  "Error accessing GDF plot.",
  "Couldnt write a file.",
  "A memory allocation failed.",
  "The anatomical volume is not loaded.",
  "Tried to call a function on an unloaded surface.",
  "Functional overlay is not loaded.",
  "GCA volume is not loaded.",
  "Segmentation data is not loaded.",
  "The FA and EV volumes are different sizes.",
  "The aux volume must be the same size as the main volume.",
  "Couldn't cache the script name.",
  "Invalid script name.",
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

static char *sPialName = "pial" ;
static char *sOrigName = "orig" ;

static int zf, ozf;
static float fsf;

static int editflag = TRUE;
static int surflinewidth = 1;
static int gbScaleUpVolume = FALSE ;
static int gbForceEnableControlPoints = FALSE;
static int gbControlPointsChanged = FALSE ;
static int gbSavingControlPoints = FALSE ;

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

/* subdirectories local to subject's home dir */
char *ksaFileNameSubDirs[tkm_knNumFileNameTypes] = {
  ".", "fmri", "mri", "bem", "surf", "mri",
  "mri/transforms", "label", "mri", "",
  "image/rgb", "tmp", "tmp", "lib/tcl", "touch"
};

/* input starts with ., gsUserHomeDir will be prepended. if ~ or nothing,
   gsSubjectHomeDir will be prepended. anything else, first use the subject
   home dir and then the subdir of the proper type, then the proper
   remainder of the input file name. if it starts with a /, nothing will be
   added. */
void MakeFileName ( char*         isInput,
                    tkm_tFileName iType,
                    char*         osCompleteFileName,
                    int           inDestSize );

/* only tries to prepend subdirectories if this flag is true. it is set to
   false if -f is used on the command line. */
tBoolean gEnableFileNameGuessing = TRUE;
tBoolean gbGuessWarningSent = FALSE;

static  tBoolean     bExit         = FALSE;


// ==========================================================================

// ================================================== SELECTING CONTROL POINTS

x3DListRef gControlPointList = NULL;

/* Conformed volumes are 256^3 so that's our size. */
#define kControlPointSpaceDimension 256

/* Inits the control point list to the dimensions of the main
   anatomical volume. */
tkm_tErr InitControlPointList ();

/* Reads the control.dat file, transforms all pts from RAS space to
   voxel space, and adds them as control pts */
void ReadControlPointFile ();

/* Writes all control points to the control.dat file in RAS space */
void WriteControlPointFile   ();

/* Returns distances to nearest control point on the same plane. returns
   0 if there isn't one. */
float FindNearestMRIIdxControlPoint ( xVoxelRef   iAnaIdx,
                                      mri_tOrientation iPlane,
                                      xVoxelRef*   opCtrlPt );

/* Makes a copy of the ctrl pt and puts it in the ctrl pt list */
void AddMRIIdxControlPoint ( xVoxelRef iMRIIdx, tBoolean ibWriteToFile );

/* Removes ctrl pt from the list and deletes it */
void DeleteMRIIdxControlPoint ( xVoxelRef iMRIIdx, tBoolean  ibWriteToFile  );

/* Function for comparing voxels */
xList_tCompare CompareVoxels ( void* inVoxelA, void* inVoxelB );

// ===========================================================================

// ================================================================== SURFACES

#include "mriSurface.h"

/* This flag toggles between using sufaceRAS coordinates and worldRAS
   coordinates for all RAS based output that is going to tksurfer at
   some point. This flag is actually specified in the surface
   strucutre, useRealRAS. If it is 1, we use worldRAS functions, if
   not, we use surfaceRAS functions. If a surface is not loaded, it's
   0 by default. If a surface is load with useRealRAS on, the user
   will be asked if they want to switch modes. */
tBoolean gbUseRealRAS = FALSE;
tBoolean gbSetFirstUseRealRAS = TRUE;

mriSurfaceRef gSurface[tkm_knNumSurfaceTypes];

tkm_tErr LoadSurface          ( tkm_tSurfaceType iType,
                                char*            isName );
tkm_tErr LoadSurfaceVertexSet ( tkm_tSurfaceType iType,
                                Surf_tVertexSet  iSet,
                                char*            fname );
void   UnloadSurface          ( tkm_tSurfaceType iType );


tkm_tErr LoadSurfaceAnnotation ( tkm_tSurfaceType iType,
                                 char*            isName );

void   WriteSurfaceValues     ( tkm_tSurfaceType iType,
                                char*            isFileName );

void GotoSurfaceVertex                    ( Surf_tVertexSet iSurface,
                                            int             inVertex );
void FindNearestSurfaceVertex             ( Surf_tVertexSet iSet );
void FindNearestInterpolatedSurfaceVertex ( Surf_tVertexSet iSet );
void AverageSurfaceVertexPositions        ( int             inNumAverages );

void SetUseRealRAS ( tkm_tSurfaceType iType,
                     tBoolean ibUseRealRAS );

/* Calculates a surface client transformation based on the
   gbUseRealRAS and sets it in the surface. */
tkm_tErr CalcAndSetSurfaceClientTransformation( tkm_tSurfaceType iType );

// ===========================================================================

// ========================================================= SELECTING REGIONS

/* selecting regions works much like editing. the user chooses a brush
   shape and size and paints in a region. the tool can be toggled
   between selecting and unselecting. the selected pixels are kept in
   a voxel space for optimized retreival in the draw loop. there are
   the usual functions for adding and removing voxels as well as
   saving them out to a file. */

mriVolumeRef gSelectionVolume = NULL;
int     gSelectionVolumeXDimension;
int     gSelectionVolumeYDimension;
int     gSelectionVolumeZDimension;
int     gSelectionCount;
xListRef gSelectionList = NULL;

tkm_tErr InitSelectionModule ();
void   DeleteSelectionModule ();

tkm_tErr AllocateSelectionVolume ();

/* adds or removes voxels to selections. if a voxel that isn't in the
   selection is told to be removed, no errors occur. this is called from the
   brush function. The iaAnaIdx parameter can be an array. */
void AddVoxelsToSelection      ( xVoxelRef  iaMRIIdx, int inCount );
void RemoveVoxelsFromSelection ( xVoxelRef  iaMRIIdx, int inCount );

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

typedef struct {
  tBoolean   mbSelect; /* 1 for select, 0 for deselect */
  int        mnCount;
}
tkm_tFloodSelectCallbackData;
tkm_tErr FloodSelect ( xVoxelRef         iSeedAnaIdx,
                       tBoolean          ib3D,
                       tkm_tVolumeTarget iSrc,
                       float             iFuzzy,
                       float             iDistance,
                       tBoolean          ibSelect );

/* Callback for the flood. */
Volm_tVisitCommand FloodSelectCallback ( xVoxelRef iMRIIdx,
                                         float     iValue,
                                         void*     iData );

// ===========================================================================

// ============================================================= VOLUME ACCESS

#include "mriVolume.h"

static mriVolumeRef     gAnatomicalVolume[tkm_knNumVolumeTypes];
static int              gnAnatomicalDimensionX = 0;
static int              gnAnatomicalDimensionY = 0;
static int              gnAnatomicalDimensionZ = 0;
static mriTransformRef  gMRIIdxToAnaIdxTransform    = NULL;
static tBoolean         gbAnatomicalVolumeDirty[tkm_knNumVolumeTypes];

tkm_tErr LoadVolume    ( tkm_tVolumeType iType,
                         char*           isFileName,
                         tBoolean        ibConform );

tkm_tErr UnloadVolume  ( tkm_tVolumeType iType );

/* saves the anatomical volume. if isPath is null, saves over the original. */
tkm_tErr SaveVolume    ( tkm_tVolumeType iType,
                         char*           isPath,
                         tBoolean        ibSetFileName );

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
                                      float           ifBrightness,
                                      float           ifContrast );
void SetVolumeColorMinMax           ( tkm_tVolumeType iVolume,
                                      float           ifMin,
                                      float           ifMax );

void SetVolumeSampleType  ( tkm_tVolumeType  iVolume,
                            Volm_tSampleType iType );

void SetVolumeResampleMethod  ( tkm_tVolumeType      iVolume,
                                Volm_tResampleMethod iMethod );

void ThresholdVolume ( int    inLevel,
                       tBoolean      ibAbove,
                       int        inNewLevel );
void FlipVolume       ( mri_tOrientation iAxis );
void RotateVolume    ( mri_tOrientation iAxis,
                       float      ifDegrees );

void EditAnatomicalVolumeInRangeArray ( tkm_tVolumeType iVolume,
                                        xVoxelRef       iaMRIIdx,
                                        int             inCount,
                                        Volm_tValue     inLow,
                                        Volm_tValue     inHigh,
                                        Volm_tValue     inNewValue );

void CloneAnatomicalVolumeInRangeArray ( tkm_tVolumeType iDestVolume,
                                         tkm_tVolumeType iSourceVolume,
                                         xVoxelRef       iaMRIIdx,
                                         int             inCount,
                                         Volm_tValue     inLow,
                                         Volm_tValue     inHigh );

void SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
                                 int             iAnaX0,
                                 int             iAnaX1,
                                 int             iAnaY0,
                                 int             iAnaY1,
                                 int             iAnaZ0,
                                 int             iAnaZ1,
                                 float           iNewValue );

int EditAnatomicalVolume ( xVoxelRef iMRIIdx, int inValue );
int EditAuxAnatomicalVolume ( xVoxelRef iMRIIdx, int inValue );

typedef struct {
  int             mnValue;
  int             mnCount;
  tkm_tVolumeType mVolume;
}
tkm_tFloodFillAnatomicalCallbackData;

tkm_tErr FloodFillAnatomicalVolume ( tkm_tVolumeType iVolume,
                                     xVoxelRef       iAnaIdx,
                                     int             inValue,
                                     tBoolean        ib3D,
                                     float           iFuzzy,
                                     float           iDistance );
/* Callback for the flood. */
Volm_tVisitCommand FloodFillAnatomicalCallback ( xVoxelRef iMRIIdx,
                                                 float     iValue,
                                                 void*     iData );
void ConvertRASToAnaIdx ( xVoxelRef iRAS,
                          xVoxelRef oAnaIdx );

void ConvertAnaIdxToRAS ( xVoxelRef iAnaIdx,
                          xVoxelRef oRAS );

float gfaBrightness[tkm_knNumVolumeTypes];
float gfaContrast[tkm_knNumVolumeTypes];
float gfaAnaColorMin[tkm_knNumVolumeTypes];
float gfaAnaColorMax[tkm_knNumVolumeTypes];

void SetVolumeBrightnessAndContrast  ( tkm_tVolumeType iVolume,
                                       float        ifBrightness,
                                       float        ifContrast );
void SendVolumeColorScaleUpdate ( tkm_tVolumeType iVolume );

void SetCursorToCenterOfVolume ( tkm_tVolumeType iVolume );

// ========================================================= FUNCTIONAL VOLUME

#include "tkmFunctionalVolume.h"

tkmFunctionalVolumeRef gFunctionalVolume = NULL;

tkm_tErr LoadFunctionalOverlay    ( char* isFileName,
                                    char* isOffsetFileName,
                                    FunD_tRegistrationType iRegType,
                                    char* isRegistrationFileName );
tkm_tErr LoadFunctionalTimeCourse ( char* isFileName,
                                    char* isOffsetFileName,
                                    FunD_tRegistrationType iRegType,
                                    char* isRegistrationFileName );

tkm_tErr SmoothOverlayData ( float ifSigma );

void TranslateOverlayRegisgtration ( float ifDistance, char isDirection );
void RotateOverlayRegistration     ( float ifDegrees, char isDirection );
void ScaleOverlayRegisgtration     ( float ifDistance, char isDirection );

// ===========================================================================

// ============================================================== SEGMENTATION

#include "colortab.h"

#define kfDefaultSegmentationAlpha 0.3

/* gSegmentationVolume is the normal segmentation volume,
   gPreviousSegmentationVolume is a backup that is created before it
   is recomputed in RecomputeSegmentation, and
   gSegmentationChangedVolume contains flags at the values that were
   changed in this session. */
static mriVolumeRef  gSegmentationVolume[tkm_knNumSegTypes];
static mriVolumeRef  gPreviousSegmentationVolume[tkm_knNumSegTypes];
static mriVolumeRef  gSegmentationChangedVolume[tkm_knNumSegTypes];

static COLOR_TABLE*  gColorTable[tkm_knNumSegTypes];
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
                                  char*           inColorFileName,
                                  tBoolean        ibConform );

/* Saves a segmentation volume to COR format. */
void SaveSegmentationVolume     ( tkm_tSegType    iVolume,
                                  char*           inFileName,
                                  tBoolean        ibSetFileName );

/* Uses the ChangedVolume to save a volume with only those label
   values that have been editied in this session. Basically an AND of
   the SegmentationVolume and ChangedVolume. */
typedef struct {
  mriVolumeRef mSrcVolume;
  mriVolumeRef mDestVolume;
}
tkm_tExportSegmentationParams, *tkm_tExportSegmentationParamsRef;
tkm_tErr ExportChangedSegmentationVolume ( tkm_tSegType    iVolume,
                                           char*           inVolumeDir );

/* Creates an segmentation from a surface annotation file. Only the
   voxel values that lie on the surface will be filled. */
tkm_tErr ImportSurfaceAnnotationToSegmentation ( tkm_tSegType iVolume,
                                                 char*    inAnnotationFileName,
                                                 char*    inColorFileName );

tkm_tErr LoadSegmentationColorTable ( tkm_tSegType iVolume,
                                      char*        inColorFileName );

/* Updates the tcl interface with a new color table. */
tkm_tErr SendColorTableInformationToTcl ( tkm_tSegType iVolume );


/* Sets the display alpha for determining the opacity of the
   segmentation overlay. */
void SetSegmentationAlpha ( float ifAlpha );

/* Gets a color at an index in a seg volume. Given a base color, this
   function will blend the label color according to the opacity alpha
   and return the new color. */
void GetSegmentationColorAtVoxel ( tkm_tSegType iVolume,
                                   xVoxelRef    iMRIIdx,
                                   xColor3fRef  iBaseColor,
                                   xColor3fRef  oColor );

/* Returns the name of the label and/or label index for a location. */
void GetSegLabel  ( tkm_tSegType iVolume,
                    xVoxelRef    iMRIIdx,
                    int*         onIndex,
                    char*        osLabel );

/* Gets the segmentation value at the cursor and calls
   SelectSegLabel. */
tkm_tErr SelectSegLabelAtCursor ( tkm_tSegType iVolume );

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
static void RecomputeUpdateCallback ( MRI* iMRIValues );

/* A call back for a visit function. If the two values are equal
   (there's a float in ipnTarget) then this voxel will be added to the
   selection. */
Volm_tVisitCommand AddSimilarVoxelToSelection ( xVoxelRef iMRIIdx,
                                                float     iValue,
                                                void*     ipnTarget );
/* A call back for a visit function. If the two values are equal
   (there's a float in ipnTarget) then this voxel will be added to the
   functional selection so it can later be graphed. */
Volm_tVisitCommand AddSimilarVoxelToGraphAvg  ( xVoxelRef iMRIIdx,
                                                float     iValue,
                                                void*     ipnTarget );

/* A call back for a visit function. When run on the changed volume by
   ExportChangedSegmentationVolume, for each value that is 1, sets the
   value in the dest volume to the value in the src volume. Uses
   tkm_tExportSegmentationParamsRef. */
Volm_tVisitCommand SetChangedSegmentationValue ( xVoxelRef iMRIIdx,
                                                 float     iValue,
                                                 void*     iData );

/* A callback to an undo entry. Sets the segmentation volume value at
   this index, as well as the changed volume index.  */
int  EditSegmentation ( tkm_tSegType iVolume,
                        xVoxelRef    iMRIIdx,
                        int          inIndex );

/* Calculates the volume of a contiguous label. */
typedef struct {
  float mSourceLabel;
  int   mCount;
}
tkm_tSumSimilarSegmentationParams;
void CalcSegLabelVolume   ( tkm_tSegType iVolume,
                            xVoxelRef    iMRIIdx,
                            int*         onVolume );
Volm_tVisitCommand SumSimilarValues ( xVoxelRef iMRIIdx,
                                      float     iValue,
                                      void*     iData );

/* Use to set the segmentation value. Also adds to the undo
   list. TODO: Add it to the undo list, it doesn't currently. */
void SetSegmentationValue    ( tkm_tSegType iVolume,
                               xVoxelRef    iMRIIdx,
                               int          inIndex );
void SetSegmentationValues   ( tkm_tSegType iVolume,
                               xVoxelRef    iaMRIIdx,
                               int          inCount,
                               int          inIndex );

/* Flood fills a segmentation volume with various parameters. */
typedef struct {
  tkm_tSegType mTargetVolume;
  int          mnNewSegLabel;
  int          mnCount;
}
tkm_tFloodFillCallbackData;
tkm_tErr FloodFillSegmentation ( tkm_tSegType      iVolume,
                                 xVoxelRef         iMRIIdx,
                                 int               inIndex,
                                 tBoolean          ib3D,
                                 tkm_tVolumeTarget iSrc,
                                 float             iFuzzy,
                                 float             iDistance );
/* Callback for the flood. */
Volm_tVisitCommand FloodFillSegmentationCallback ( xVoxelRef iMRIIdx,
                                                   float     iValue,
                                                   void*     iData );

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
}
UndoEntry, *UndoEntryRef;

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

int            gUndoVolumeDimensions[3];
tBoolean***    gUndoVoxelFlagVolume = NULL;
Volm_tValue*** gUndoVoxelValueVolume = NULL;

typedef struct {
  int mRestoreValue;
}
UndoVolumeEntry, *UndoVolumeEntryRef;

/* initializes and deletes the volumes */
tkm_tErr InitUndoVolume    ();
void   DeleteUndoVolume ();

/* adds a value to the volume at an anatomical index */
void   AddMRIIdxAndValueToUndoVolume ( xVoxelRef    iMRIIdx,
                                       Volm_tValue  iValue );

/* sees if there is a value for this MRI idx, i.e. if it can be undone */
tBoolean IsMRIIdxInUndoVolume         ( xVoxelRef iMRIIdx );

/* resotres the values for all voxels touching this one */
void   RestoreUndoVolumeAroundMRIIdx ( xVoxelRef iMRIIdx );

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

void AlignSelectedHeadPointToMRIIdx ( xVoxelRef iMRIIdx );

/* ======================================================================= */

/* =================================================================== GCA */

#include "gca.h"

GCA* gGCAVolume     = NULL;
MRI *gMRI = NULL ;
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

/* ============================================================ GDF VOLUME */

int gGDFID = -1;
mriTransformRef gMRIIdxToGDFIdxTransform = NULL;

tkm_tErr LoadGDFHeader ( char* isFSGDHeader,
                         int iRegistrationType,
                         char* isRegistrationFile );

tkm_tErr GDFPlotMRIIdx ( xVoxelRef iMRIIdx );

tkm_tErr BeginGDFPointList ();
tkm_tErr AddGDFPlotMRIIdx ( xVoxelRef iMRIIdx );
tkm_tErr EndGDFPointList ();

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

/* Find a good edit.dat file. */
void CopyEditDatFileName ( char* osFileName, int izFileName );

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
char* SendTCLCommand ( char * inCommand );


xGrowableArrayRef gCachedTclCommands = NULL;
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

int LoadOrigSurf = 1;
int LoadPialSurf = 1;

// ===========================================================================

#ifdef Linux
extern void scale2x(int, int, unsigned char *);
#endif

void ParseCmdLineArgs( int argc, char *argv[] );
void WriteVoxelToEditFile    ( xVoxelRef inVolumeVox );

void rotate_brain(float a,char c) ;
void translate_brain(float a,char c) ;
void UpdateAndRedraw ();
void pix_to_rgb(char *fname) ; // another method of saving a screenshot
void scrsave_to_rgb(char *fname) ; // use scrsave to save a screen shot
void save_rgb(char *fname) ; // saves a screen shot
void save_tiff(char *fname) ;
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
void checkLicense(char* dirname) {
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
  if (lfile) {
    fscanf(lfile,"%s\n",email);
    fscanf(lfile,"%s\n",magic);
    fscanf(lfile,"%s\n",key);

    sprintf(gkey,"%s.%s",email,magic);
    if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0) {
      printf("No valid license key !\n");
      exit(-1);
    }
  } else {
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

  tkm_tErr     eResult                    = tkm_tErr_NoErr;
  FunV_tErr    eFunctional                = FunV_tErr_NoError;
  int          nNumProcessedVersionArgs   = 0;
  int          nCurrentArg                = 0;
  int          bFatalError                = FALSE;
  char         sArg[tkm_knPathLen]        = "";
  char         sError[tkm_knErrStringLen] = "";
  char*        pEnvVar                    = NULL;

  tBoolean     bSubjectDeclared   = FALSE;
  tBoolean     bGetSubjectFromReg = FALSE;
  tBoolean     bUsingMRIRead      = FALSE;

  char         sSubject[tkm_knPathLen]  = "";
  char         sImageDir[tkm_knPathLen] = "";
  tBoolean     bScaleUpVolume           = FALSE;
  tBoolean     bConformVolume           = FALSE;

  tBoolean     bSurfaceDeclared        = FALSE;
  char         sSurface[tkm_knPathLen] = "";

  tBoolean     bAuxSurfaceDeclared       = FALSE;
  char         sAuxSurface[tkm_knPathLen] = "";

  tBoolean     bLocalImageDir = FALSE;
  tBoolean     bNoEdit        = FALSE;

  tBoolean     bLoadingAuxVolume         = FALSE;
  char         sAuxVolume[tkm_knPathLen] = "";
  tBoolean     bConformAuxVolume         = FALSE;

  tBoolean     bLoadingMainTransform         = FALSE;
  char         sMainTransform[tkm_knPathLen] = "";
  tBoolean     bLoadingAuxTransform          = FALSE;
  char         sAuxTransform[tkm_knPathLen]  = "";

  tBoolean               bLoadingOverlay                      = FALSE;
  char                   sOverlayFileName[tkm_knPathLen]      = "";
  char                   sOverlayOffsetFileName[tkm_knPathLen]= "";
  FunD_tRegistrationType overlayRegType             = FunD_tRegistration_None;
  char                   sOverlayRegistration[tkm_knPathLen]  = "";
  char*                  psOverlayOffsetFileName              = NULL;

  tBoolean               bLoadingTimeCourse                       = FALSE;
  char                   sTimeCourseFileName[tkm_knPathLen]       = "";
  char*                  psTimeCourseOffsetFileName               = NULL;
  char                   sTimeCourseOffsetFileName[tkm_knPathLen] = "";
  FunD_tRegistrationType timecourseRegType          = FunD_tRegistration_None;
  char                   sTimeCourseRegistration[tkm_knPathLen]   = "";

  tBoolean     bEnablingRegistration      = FALSE;

  tBoolean     bLoadingSegmentation                  = FALSE;
  char         sSegmentationPath[tkm_knPathLen]      = "";
  char         sSegmentationColorFile[tkm_knPathLen] = "";
  tBoolean     bConformSegmentation                  = FALSE;

  tBoolean     bLoadingAuxSegmentation                  = FALSE;
  char         sAuxSegmentationPath[tkm_knPathLen]      = "";
  char         sAuxSegmentationColorFile[tkm_knPathLen] = "";
  tBoolean     bConformAuxSegmentation                  = FALSE;

  tBoolean     bLoadingVLI              = FALSE;
  char         sVLIFile1[tkm_knPathLen] = "";
  char         sVLIFile2[tkm_knPathLen] = "";

  tBoolean              bThresh              = FALSE;
  FunV_tFunctionalValue min                  = 0;
  tBoolean              bMid                 = FALSE;
  FunV_tFunctionalValue mid                  = 0;
  tBoolean              bSlope               = FALSE;
  FunV_tFunctionalValue slope                = 0;
  tBoolean              bMax                 = FALSE;
  FunV_tFunctionalValue max                  = 0;

  tBoolean              bSmooth              = FALSE;
  float                 smoothSigma          = 0;
  tBoolean              bRevPhaseFlag        = FALSE;
  int                   nRevPhaseFlag        = 0;
  tBoolean              bTruncPhaseFlag      = FALSE;
  int                   nTruncPhaseFlag      = 0;
  tBoolean              bUseOverlayCacheFlag = FALSE;
  int                   nUseOverlayCacheFlag = 0;

  tBoolean      bLoadingHeadPts                  = FALSE;
  tBoolean      bHaveHeadPtsTransform            = FALSE;
  char          sHeadPts[tkm_knPathLen]          = "";
  char          sHeadPtsTransform[tkm_knPathLen] = "";

  char          sLabel[tkm_knPathLen] = "";
  tBoolean      bLoadingLabel         = FALSE;

  char          sSubjectTest[tkm_knPathLen] = "";
  char          sScriptName[tkm_knPathLen]  = "";

  float         fSegmentationAlpha      = 0;
  tBoolean      bSegmentationAlpha      = FALSE;

  float         fBrightnessMain       = 0;
  float         fContrastMain         = 0;
  tBoolean      bBrightContrastMain   = FALSE;
  float         fBrightnessAux        = 0;
  float         fContrastAux          = 0;
  tBoolean      bBrightContrastAux    = FALSE;

  float         fColorMinMain  = 0;
  float         fColorMaxMain  = 0;
  tBoolean      bColorMain     = FALSE;
  float         fColorMinAux   = 0;
  float         fColorMaxAux   = 0;
  tBoolean      bColorAux      = FALSE;

  tBoolean      bLoadingAnnotation = FALSE;
  char          sAnnotation[tkm_knPathLen] = "";
  char          sAnnotationColorTable[tkm_knPathLen] = "";

  tBoolean      bMIP = FALSE;
  char tmpstr[2000];

  DebugEnterFunction( ("ParseCmdLineArgs( argc=%d, argv=%s )",
                       argc, argv[0]) );

  /* first get the functional threshold so we don't overwrite the defaults */
  DebugNote( ("Getting default functional threshold") );
  eFunctional = FunV_GetThreshold( gFunctionalVolume, &min, &mid, &slope );
  DebugAssertThrow( (FunV_tErr_NoError == eFunctional ) );

  // Read in env defaults
  if(getenv("FS_TKFTHRESH")){
    sscanf(getenv("FS_TKFTHRESH"),"%f",&min);
    bThresh = TRUE;
  }
  if(getenv("FS_TKFMAX")){
    sscanf(getenv("FS_TKFMAX"),"%f",&max);
    bMax = TRUE;
  }

  /* First look for the version option and handle that. If found,
     shorten our argc and argv count. If those are the only args we
     had, exit. */
  /* rkt: check for and handle version tag */
  nNumProcessedVersionArgs =
    handle_version_option
    (argc, argv,
     "$Id: tkmedit.c,v 1.346 2011/10/15 15:43:39 fischl Exp $",
     "$Name:  $");
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
    printf("subject image_type  : looks in $SUBJECTS_DIR/subject/mri "
           "for image_type. If subject=getreg, gets subject from ov reg\n");
    printf("-f absolute_path    : specify volume directory or file\n");
    printf("\n");
    printf("Surface\n");
    printf("\n");
    printf("surface_file   : surface file to load (relative "
           "to $SUBJECTS_DIR/surf \n");
    printf("               : or absolute)\n");
    printf("\n");
    printf("Options\n");
    printf("\n");
    printf("-aux <volume>  : load volume as auxilliary anatomical volume. "
           "relative to\n");
    printf("               : in $SUBJECTS_DIR/subject/mri or "
           "specify absolute path\n");
    printf("\n");
    printf("-pial <surface>        : load pial surface locations from  <surface>\n") ;
    printf("-orig <surface>        : load orig surface locations from  <surface>\n") ;
    printf("-surface <surface>     : load surface as main surface. Relative to\n");
    printf("                       : in $SUBJECTS_DIR/subject/surf "
           "or specify absolute path\n");
    printf("-aux-surface <surface> : load surface as auxilliary surface. "
           "relative to\n");
    printf("                       : in $SUBJECTS_DIR/subject/surf or "
           "specify absolute path\n");
    printf("-surfs : load lh.white and rh.white\n");
    printf("-annotation <annotation> <color table> : import a "
           "surface annotation and color\n");
    printf("                                         table to "
           "the main surface\n");
    printf("\n");
    printf("-main-transform <transform> : loads a display transform "
           "for the main volume\n");
    printf("-aux-transform <transform>  : loads a display transform "
           "for the aux volume\n");
    printf("\n");
    printf("-bc-main <brightness> <contrast> : brightness and contrast "
           "for main volume\n");
    printf("-bc-main-fsavg use .58 and 14 for main volume (good for fsaverage)\n");
    printf("-mm-main <min> <max>             : color scale min and "
           "max for main volume\n");
    printf("-bc-aux <brightness> <contrast>  : brightness and "
           "contrast for aux volume\n");
    printf("-mm-aux <min> <max>              : color scale min "
           "and max for aux volume\n");
    printf("\n");
    printf("-conform      : conform the main anatomical volume\n");
    printf("-aux-conform  : conform the aux anatomical volume\n");
    printf("\n");
    printf("-overlay <file>            : load functional overlay volume (-ov)\n");
    printf("-overlay-reg-find          : find overlay registration (-ovreg)"
           "volume in data dir\n");
    printf("-overlay-reg-identity      : generate identity for "
           "overlay registration volume\n");
    printf("-overlay-reg <registration> : load registration file "
           "for overlay volume \n");
    printf("                            : (default is register.dat "
           "in same path as\n");
    printf("                            :  volume)\n");
    printf("\n");
    printf("-reg regfile : use regfile for both overlay and time course\n");
    printf("\n");
    printf("-fthresh <value> : threshold for overlay (FS_TKFTHRESH)\n");
    printf("-fmax <value>    : max/sat for overlay (FS_TKFMAX)\n");
    printf("-fmid <value>          : values for functional overlay display\n");
    printf("-fslope <value>        : (default is 0, 1.0, and 1.0)\n");
    printf("-fsmooth <sigma>       : smooth functional overlay "
           "after loading\n");
    printf("-noblend               : do not blend activation color with backgound\n");
    printf("\n");
    printf("-revphaseflag <1|0>      : reverses phase display"
           " in overlay (default off)\n");
    printf("-truncphaseflag <1|0>    : truncates overlay "
           "values below 0 (default off)\n");
    printf("-overlaycache <1|0>      : uses overlay cache (default off)\n");
    printf("\n");
    printf("-sdir <subjects dir>       : (default is getenv(SUBJECTS_DIR)\n");
    printf("-timecourse <file>         : load functional timecourse volume (-t)\n");
    printf("-timecourse-reg-find       : find timecourse registration (-treg)"
           "volume in data dir\n");
    printf("-timecourse-reg-identity   : generate identity for "
           "timecourse registration volume\n");
    printf("-timecourse-reg <registration> : load registration "
           "file for timecourse   \n");
    printf("                               : volume (default "
           "is register.dat in\n");
    printf("                               : same path as volume)\n");
    printf("-timecourse-offset <path/stem> : load timecourse offset volume\n");
    printf("\n");
    printf("-segmentation <volume> [colors]    : load segmentation "
           "volume and color file\n");
    printf("-aux-segmentation <volume> [colors]: load aux "
           "segmentation volume and color file\n");
    printf("-segmentation-opacity <opacity>    : opacity of the "
           "segmentation \n");
    printf("                                   : overlay (default is 0.3)\n");
    printf("-aseg : load aseg.mgz and standard color table\n");
    printf("-save-cpts : automatically save control point file after adding or deleting\n") ;
    printf("-aparc+aseg : load aparc+aseg.mgz and standard color table\n");
    printf("-wmparc : load wmparc.mgz and standard color table\n");
    printf("\n");
    printf("-seg-conform      : conform the main segmentation volume\n");
    printf("-aux-seg-conform  : conform the aux segmentation volume\n");
    printf("\n");
    printf("-headpts <points> [<trans>]   : load head points file "
           "and optional\n");
    printf("                              : transformation\n");
    printf("\n");
    printf("-mip: turn on maximum intensity projection");
    printf("\n");
    printf("-interface script    : specify interface script "
           "(default is tkmedit.tcl)\n");
    printf("\n");
    printf("-exit                : exit after rendering\n");
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
  while ( nCurrentArg < argc
          && FALSE == bFatalError) {

    /* get the arg */
    DebugNote( ("Copying argv[%d]", nCurrentArg) );
    xUtil_strncpy( sArg, argv[nCurrentArg], sizeof(sArg) );

    /* check for a option */
    if ( '-' == sArg[0] ) {

      if ( MATCH( sArg, "-tcl" ) ) {

        /* check for one more. */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-bc-main" ) ) {

        /* check for the 2 values following the switch */
        if ( argc > nCurrentArg + 2
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

      } else if (MATCH(sArg, "-bc-main-fsavg"  ) ) {
        DebugNote( ("Parsing -bc-main-fsavg option") );
        fBrightnessMain = .58;
        fContrastMain = 14;
        bBrightContrastMain = TRUE;
        nCurrentArg ++ ;

      } else if (MATCH(sArg, "-scaleup"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling scaleup.") );
        gbScaleUpVolume = bScaleUpVolume = TRUE;
        nCurrentArg ++;

      } else if (MATCH(sArg, "-enablecontrolpoints"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling control points..") );
        gbForceEnableControlPoints = TRUE;
        nCurrentArg ++;

      } else if (MATCH(sArg, "-conform"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling conform.") );
        bConformVolume = TRUE;
        nCurrentArg ++;

      } else if ( MATCH( sArg, "-mm-main" ) ) {

        /* check for the 2 values following the switch */
        if ( argc > nCurrentArg + 2
             && '-' != argv[nCurrentArg+1][0]
             && '-' != argv[nCurrentArg+2][0] ) {

          /* get the values */
          DebugNote( ("Parsing -mm-main option") );
          fColorMinMain = atof( argv[nCurrentArg+1] );
          fColorMaxMain = atof( argv[nCurrentArg+2] );
          bColorMain = TRUE;
          nCurrentArg +=3 ;

        } else {

          /* misuse of that switch */

          tkm_DisplayError( "Parsing -mm-main option",
                            "Expected two arguments",
                            "This option needs two arguments: the min"
                            "and the max color values for the volume." );
          nCurrentArg += 1;
        }

      } else if ( MATCH( sArg, "-bc-aux" ) ) {

        /* check for the 2 values following the switch */
        if ( argc > nCurrentArg + 2
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

      } else if ( MATCH( sArg, "-mm-aux" ) ) {

        /* check for the 2 values following the switch */
        if ( argc > nCurrentArg + 2
             && '-' != argv[nCurrentArg+1][0]
             && '-' != argv[nCurrentArg+2][0] ) {

          /* get the values */
          DebugNote( ("Parsing -mm-aux option") );
          fColorMinAux = atof( argv[nCurrentArg+1] );
          fColorMaxAux = atof( argv[nCurrentArg+2] );
          bColorAux = TRUE;
          nCurrentArg +=3 ;

        } else {

          /* misuse of that switch */
          tkm_DisplayError( "Parsing -mm-aux option",
                            "Expected two arguments",
                            "This option needs two arguments: the min"
                            "and the max color values for the volume." );
          nCurrentArg += 1;
        }

      } 
      else if ( MATCH( sArg, "-surfs" ) ) {
	xUtil_strncpy( sSurface, "lh.white", sizeof(sSurface) );
	bSurfaceDeclared = TRUE;
	xUtil_strncpy( sAuxSurface, "rh.white", sizeof(sAuxSurface) );
	bAuxSurfaceDeclared = TRUE;
	LoadOrigSurf = 0;
	nCurrentArg ++;
      }
      else if ( MATCH( sArg, "-defects" ) || MATCH( sArg, "-lh-defects" ) ) {
	xUtil_strncpy( sSurface, "lh.orig", sizeof(sSurface) );//green
	bSurfaceDeclared = TRUE;
	xUtil_strncpy( sAuxSurface, "lh.orig.nofix", sizeof(sAuxSurface) );//yellow
	bAuxSurfaceDeclared = TRUE;
	LoadOrigSurf = 1;
	LoadPialSurf = 0;
        xUtil_strncpy( sSegmentationPath, "surface.defects.mgz",
                       sizeof(sSegmentationPath) );
        pEnvVar = getenv("FREESURFER_HOME");
        sprintf( sSegmentationColorFile,"%s/DefectLUT.txt", pEnvVar );
        bLoadingSegmentation = TRUE;
	fSegmentationAlpha = 0.8;
	bSegmentationAlpha = TRUE;
	xUtil_strncpy( sAuxVolume,"wm.mgz",sizeof(sAuxVolume) );
	bLoadingAuxVolume = TRUE;
	nCurrentArg ++;
      }
      else if ( MATCH( sArg, "-rh-defects" ) ) {
	xUtil_strncpy( sSurface, "rh.orig", sizeof(sSurface) );//green
	bSurfaceDeclared = TRUE;
	xUtil_strncpy( sAuxSurface, "rh.orig.nofix", sizeof(sAuxSurface) );//yellow
	bAuxSurfaceDeclared = TRUE;
	LoadOrigSurf = 1;
	LoadPialSurf = 0;
        xUtil_strncpy( sSegmentationPath, "surface.defects.mgz",
                       sizeof(sSegmentationPath) );
        pEnvVar = getenv("FREESURFER_HOME");
        sprintf( sSegmentationColorFile,"%s/DefectLUT.txt", pEnvVar );
        bLoadingSegmentation = TRUE;
	fSegmentationAlpha = 0.8;
	bSegmentationAlpha = TRUE;
	xUtil_strncpy( sAuxVolume,"wm.mgz",sizeof(sAuxVolume) );
	bLoadingAuxVolume = TRUE;
	nCurrentArg ++;
      }
      else if ( MATCH( sArg, "-aux-surface" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg into a destructible string */
          DebugNote( ("Parsing -aux-surface name") );
          xUtil_strncpy( sAuxSurface, argv[nCurrentArg+1],
                         sizeof(sAuxSurface) );
          bAuxSurfaceDeclared = TRUE;
          nCurrentArg += 2;

        } 
	else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -aux-surface option",
                            "Expected an argument",
                            "This option needs an argument, the file name "
                            "of the surface to load." );
          nCurrentArg ++;
        }

      } else if ( MATCH( sArg, "-surface" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg into a destructible string */
          DebugNote( ("Parsing -surface name") );
          xUtil_strncpy( sSurface, argv[nCurrentArg+1],
                         sizeof(sSurface) );
          bSurfaceDeclared = TRUE;
          nCurrentArg += 2;

        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -surface option",
                            "Expected an argument",
                            "This option needs an argument, the file name "
                            "of the surface to load." );
          nCurrentArg ++;
        }

      } else if ( MATCH( sArg, "-annotation" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] &&
             '-' != argv[nCurrentArg+2][0] ) {

          /* copy args */
          DebugNote( ("Parsing -annotation name") );
          xUtil_strncpy( sAnnotation, argv[nCurrentArg+1],
                         sizeof(sAnnotation) );
          xUtil_strncpy( sAnnotationColorTable, argv[nCurrentArg+2],
                         sizeof(sAnnotationColorTable) );
          bLoadingAnnotation = TRUE;
          nCurrentArg += 3;

        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -annotation option",
                            "Expected two arguments",
                            "This option needs two arguments, the file name "
                            "of the surface to load and the file name "
                            "of the color table to use." );
          nCurrentArg ++;
        }

      } else if ( MATCH( sArg, "-o" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 2 &&
             '-' != argv[nCurrentArg+1][0] &&
             '-' != argv[nCurrentArg+2][0] ) {

          /* read the overlay path and stem */
          DebugNote( ("Parsing overlay in -o option") );
          xUtil_snprintf( sOverlayFileName, sizeof(sOverlayFileName),
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
        if ( argc > nCurrentArg + 1
             && '-' != argv[nCurrentArg][0] ) {

          /* read in time course path and stem. */
          DebugNote( ("Parsing time course in -o option") );
          xUtil_snprintf( sTimeCourseFileName,
                          sizeof(sTimeCourseFileName),
                          "%s/%s",
                          argv[nCurrentArg], argv[nCurrentArg+1] );
          if (!fio_FileExistsReadable(sTimeCourseFileName)) {
            printf("ERROR: cannot find time course %s\n",sTimeCourseFileName);
            exit(1);
          }
          bLoadingTimeCourse = TRUE;
          nCurrentArg += 2;
        }

      } else if ( MATCH( sArg, "-overlay" ) || MATCH( sArg, "-ov" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg into a destructible string */
          DebugNote( ("Parsing -overlay option") );
          xUtil_strncpy( sOverlayFileName, argv[nCurrentArg+1],
                         sizeof(sOverlayFileName) );
          if (!fio_FileExistsReadable(sOverlayFileName)) {
            printf("ERROR: cannot find overlay %s\n",sOverlayFileName);
            exit(1);
          }
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

      } 
      else if ( MATCH( sArg, "-overlay-reg" ) || MATCH( sArg, "-orf" ) ||
                  MATCH( sArg, "-oreg" ) || MATCH( sArg, "-ovreg" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg  */
          DebugNote( ("Parsing -overlay-reg option") );
          xUtil_strncpy( sOverlayRegistration, argv[nCurrentArg+1],
                         sizeof(sOverlayRegistration) );
          if (!fio_FileExistsReadable(sOverlayRegistration)) {
            printf("ERROR: cannot find overlay reg %s\n",sOverlayRegistration);
            exit(1);
          }
          overlayRegType = FunD_tRegistration_File;
          nCurrentArg += 2;
          if(bGetSubjectFromReg){
            FILE *fp;
            fp = fopen(sOverlayRegistration,"r");
            fscanf(fp,"%s",sSubject);
            fclose(fp);
            //            printf("Setting subject name to %s\n",sSubject);
            DebugNote( ("Setting subject home from env") );
            eResult = SetSubjectHomeDirFromEnv( sSubject );
            DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
          }
        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -overlay-reg option",
                            "Expected an argument",
                            "This option needs an argument, the file name "
                            "of the registration data." );
          nCurrentArg ++;
        }

      } 

      //---------------------------------------------------------------
      else if ( MATCH( sArg, "-reg" ) || 
		MATCH( sArg, "-ovtreg" ) || MATCH( sArg, "-tovreg" )){
        /* make sure there are enough args */
        if( !( argc > nCurrentArg + 1 && '-' != argv[nCurrentArg+1][0]) ){
          printf("ERROR: Parsing -ovtreg option. Expected an argument.\n");
	  exit(1);
	}
	DebugNote( ("Parsing -reg option") );
	xUtil_strncpy( sOverlayRegistration, argv[nCurrentArg+1],
		       sizeof(sOverlayRegistration) );
	if (!fio_FileExistsReadable(sOverlayRegistration)) {
	  printf("ERROR: cannot find overlay reg %s\n",sOverlayRegistration);
	  exit(1);
	}
	overlayRegType = FunD_tRegistration_File;
	xUtil_strncpy( sTimeCourseRegistration, argv[nCurrentArg+1],
		       sizeof(sTimeCourseRegistration) );
	timecourseRegType = FunD_tRegistration_File;
	nCurrentArg += 2;
	if(bGetSubjectFromReg){
	  FILE *fp;
	  fp = fopen(sOverlayRegistration,"r");
	  fscanf(fp,"%s",sSubject);
	  fclose(fp);
	  printf("Setting subject name to %s\n",sSubject);
	  DebugNote( ("Setting subject home from env") );
	  eResult = SetSubjectHomeDirFromEnv( sSubject );
	  DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
	}
      } 

      else if ( MATCH( sArg, "-overlay-reg-find" ) ) {

        /* set our reg type  */
        overlayRegType = FunD_tRegistration_Find;
        nCurrentArg += 1;

      } else if ( MATCH( sArg, "-overlay-reg-identity" ) ) {

        /* set our reg type */
        overlayRegType = FunD_tRegistration_Identity;
        nCurrentArg += 1;

      } else if ( MATCH( sArg, "-overlay-offset" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg  */
          DebugNote( ("Parsing -overlay-offset option") );
          psOverlayOffsetFileName = sOverlayOffsetFileName;
          xUtil_strncpy( sOverlayOffsetFileName, argv[nCurrentArg+1],
                         sizeof(sOverlayOffsetFileName) );
          nCurrentArg += 2;

        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -overlay-offset option",
                            "Expected an argument",
                            "This option needs an argument, the path "
                            "and stem of the offset volume." );
          nCurrentArg ++;
        }

      } else if ( MATCH( sArg, "-timecourse" ) || MATCH( sArg, "-t" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 ) {

          /* copy arg into a destructible string */
          DebugNote( ("Parsing -timecourse option") );
          xUtil_strncpy( sTimeCourseFileName, argv[nCurrentArg+1],
                         sizeof(sTimeCourseFileName) );
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

      } else if ( MATCH( sArg, "-timecourse-reg" ) || 
                  MATCH( sArg, "-treg" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg  */
          DebugNote( ("Parsing -timecourse-reg option") );
          xUtil_strncpy( sTimeCourseRegistration, argv[nCurrentArg+1],
                         sizeof(sTimeCourseRegistration) );
          if (!fio_FileExistsReadable(sTimeCourseRegistration)) {
            printf("ERROR: cannot find timecourse reg %s\n",
                   sTimeCourseRegistration);
            exit(1);
          }
          timecourseRegType = FunD_tRegistration_File;
          nCurrentArg += 2;

        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -timecourse-reg option",
                            "Expected an argument",
                            "This option needs an argument, the file name "
                            "of the registration data." );
          nCurrentArg ++;
        }

      } else if ( MATCH( sArg, "-timecourse-reg-find" ) ) {

        /* set our reg type  */
        timecourseRegType = FunD_tRegistration_Find;
        nCurrentArg += 1;

      } else if ( MATCH( sArg, "-timecourse-reg-identity" ) ) {

        /* set our reg type  */
        timecourseRegType = FunD_tRegistration_Identity;
        nCurrentArg += 1;

      } else if ( MATCH( sArg, "-timecourse-offset" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy arg  */
          DebugNote( ("Parsing -timecourse-offset option") );
          psTimeCourseOffsetFileName = sTimeCourseOffsetFileName;
          xUtil_strncpy( sTimeCourseOffsetFileName, argv[nCurrentArg+1],
                         sizeof(sTimeCourseOffsetFileName) );
          nCurrentArg += 2;

        } else {

          /* misuse of that option */
          tkm_DisplayError( "Parsing -timecourse-offset option",
                            "Expected an argument",
                            "This option needs an argument, the path "
                            "and stem to the offset volume." );
          nCurrentArg ++;
        }

      } 
      else if ( MATCH( sArg, "-mni152reg" ) ){
	sprintf(sOverlayRegistration,"%s/average/mni152.register.dat",
		getenv("FREESURFER_HOME"));
	sprintf(sTimeCourseRegistration,"%s/average/mni152.register.dat",
		getenv("FREESURFER_HOME"));
	overlayRegType = FunD_tRegistration_File;
	timecourseRegType = FunD_tRegistration_File;
	nCurrentArg += 1;
      } 
      else if ( MATCH( sArg, "-register" ) ) {

        /* set our flag */
        DebugNote( ("Enabling registration.") );
        bEnablingRegistration = TRUE;
        nCurrentArg ++;

      } 
      else if ( MATCH( sArg, "-aseg" ) ) {
        xUtil_strncpy( sSegmentationPath, "aseg.mgz",
                       sizeof(sSegmentationPath) );
        pEnvVar = getenv("FREESURFER_HOME");
        sprintf( sSegmentationColorFile,"%s/FreeSurferColorLUT.txt", pEnvVar );
        bLoadingSegmentation = TRUE;
        nCurrentArg += 1;

      } 
      else if ( MATCH( sArg, "-save-cpts" ) ) {
	gbSavingControlPoints = TRUE ;
        nCurrentArg += 1;
      } 
      else if ( MATCH( sArg, "-aparc+aseg" ) ) {
        xUtil_strncpy( sSegmentationPath, "aparc+aseg.mgz",
                       sizeof(sSegmentationPath) );
        pEnvVar = getenv("FREESURFER_HOME");
        sprintf( sSegmentationColorFile,"%s/FreeSurferColorLUT.txt", pEnvVar );
        bLoadingSegmentation = TRUE;
        nCurrentArg += 1;
      } 
      else if ( MATCH( sArg, "-wmparc" ) ) {
        xUtil_strncpy( sSegmentationPath, "wmparc.mgz",
                       sizeof(sSegmentationPath) );
        pEnvVar = getenv("FREESURFER_HOME");
        sprintf( sSegmentationColorFile,"%s/FreeSurferColorLUT.txt", pEnvVar );
        bLoadingSegmentation = TRUE;
        nCurrentArg += 1;
      } 
      else if ( MATCH( sArg, "-segmentation" ) ||
                  MATCH( sArg, "-seg" ) ||
                  MATCH( sArg, "-parcellation" ) ||
                  MATCH( sArg, "-parc" ) ) {

        /* make sure there is at least one argument. */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy the filename */
          DebugNote( ("Parsing -segmentation option") );
          xUtil_strncpy( sSegmentationPath, argv[nCurrentArg+1],
                         sizeof(sSegmentationPath) );
          nCurrentArg += 2;

          /* Check for the additional LUT filename. */
          if ( argc > nCurrentArg &&
               '-' != argv[nCurrentArg][0] ) {

            /* Get the LUT from the command line. */
            xUtil_strncpy( sSegmentationColorFile, argv[nCurrentArg],
                           sizeof(sSegmentationColorFile) );
            nCurrentArg += 1;

          } else {

            /* If not, use a default color table name. */
            pEnvVar = getenv("FREESURFER_HOME");
            if ( NULL != pEnvVar ) {
              sprintf( sSegmentationColorFile,
                       "%s/FreeSurferColorLUT.txt", pEnvVar );
            } else {
              xUtil_strncpy( sSegmentationColorFile,"./FreeSurferColorLUT.txt",
                             sizeof(sSegmentationColorFile) );
            }
          }

          bLoadingSegmentation = TRUE;

        } else {

          /* misuse of that switch */
          tkm_DisplayError( "Parsing -segmentation/seg option",
                            "Expected one or two arguments",
                            "This option needs at least one argument: the "
                            "file name of the segmentation volume, and "
                            "an optional color table file name." );
          nCurrentArg ++;
        }

      } else if (MATCH(sArg, "-seg-conform"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling seg conform.") );
        bConformSegmentation = TRUE;
        nCurrentArg ++;

      } else if ( MATCH( sArg, "-aux-segmentation" ) ||
                  MATCH( sArg, "-aux-seg" ) ) {

        /* make sure there is at least one arg */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* copy the filename */
          DebugNote( ("Parsing -aux-segmentation option") );
          xUtil_strncpy( sAuxSegmentationPath, argv[nCurrentArg+1],
                         sizeof(sAuxSegmentationPath) );
          nCurrentArg += 2;

          /* Check for the additional LUT filename. */
          if ( argc > nCurrentArg &&
               '-' != argv[nCurrentArg][0] ) {

            /* Get the LUT from the command line. */
            xUtil_strncpy( sAuxSegmentationColorFile, argv[nCurrentArg],
                           sizeof(sAuxSegmentationColorFile) );
            nCurrentArg += 1;

          } else {

            /* If not, use a default color table name. */
            pEnvVar = getenv("FREESURFER_HOME");
            if ( NULL != pEnvVar ) {
              sprintf( sAuxSegmentationColorFile,
                       "%s/FreeSurferColorLUT.txt", pEnvVar );
            } else {
              xUtil_strncpy( sAuxSegmentationColorFile,
                             "./FreeSurferColorLUT.txt",
                             sizeof(sAuxSegmentationColorFile) );
            }
          }

          bLoadingAuxSegmentation = TRUE;

        } else {

          /* misuse of that switch */
          tkm_DisplayError( "Parsing -aux-segmentation/aux-seg option",
                            "Expected one or two arguments",
                            "This option needs at least one argument: the "
                            "file name of the segmentation volume, and "
                            "an optional color table file name." );
          nCurrentArg ++;
        }

      } else if (MATCH(sArg, "-aux-seg-conform"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling conform.") );
        bConformAuxSegmentation = TRUE;
        nCurrentArg ++;

      } else if ( MATCH( sArg, "-voxel-label" ) ) {

        /* make sure there are enough args */
        if ( argc > nCurrentArg + 2 &&
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

      } else if ( MATCH( sArg, "-f" ) ) {

        /* make sure subject is not already declared */
        if ( bSubjectDeclared ) {
          tkm_DisplayError( "Parsing -f option",
                            "Subject already declared",
                            "The -f option is only to be used if the "
                            "subject is not to be declared using the "
                            "subect/image type format." );
          bFatalError = TRUE;

          /* check for path */
        } else if ( argc > nCurrentArg + 1 &&
                    '-' != argv[nCurrentArg+1][0] ) {

          /* read the path */
          DebugNote( ("Parsing -f option") );
          xUtil_strncpy( sSubject, argv[nCurrentArg+1], sizeof(sSubject) );
          bUsingMRIRead = TRUE;
          bSubjectDeclared = TRUE;
          nCurrentArg += 2;

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

      } 
      else if ( MATCH( sArg, "-fminmax" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 2 &&
             '-' != argv[nCurrentArg+1][0] ) {
          /* get the value */
          DebugNote( ("Parsing -fminmax option") );
          min = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
          bThresh = TRUE;
          max = (FunV_tFunctionalValue) atof( argv[nCurrentArg+2] );
          bMax = TRUE;
          nCurrentArg +=3 ;
        } else {
          /* misuse of that switch */
          tkm_DisplayError( "Parsing -fminmax option",
                            "Expected two arguments",
                            "This option needs two arguments: the threshold and max");
          nCurrentArg += 1;
        }

      } 
      else if ( MATCH( sArg, "-fthresh" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
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

      } 
      else if ( MATCH( sArg, "-fmax" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {
          /* get the value */
          DebugNote( ("Parsing -fmax option") );
          max = (FunV_tFunctionalValue) atof( argv[nCurrentArg+1] );
          bMax = TRUE;
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

      } else if ( MATCH( sArg, "-fmid" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
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

      } else if ( MATCH( sArg, "-fslope" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
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

      } else if ( MATCH( sArg, "-sdir" ) ) {
        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
             '-' != argv[nCurrentArg+1][0] ) {

          /* get the value */
          DebugNote( ("Parsing -sdir option") );
          gsCommandLineSubjectsDir = argv[nCurrentArg+1] ;
          nCurrentArg +=2 ;

        } else {

          /* misuse of that switch */
          tkm_DisplayError( "Parsing -sdir option",
                            "Expected an argument",
                            "This option needs an argument: the subjects "
                            "dir to use." );
          nCurrentArg += 1;
        }
      } else if ( MATCH( sArg, "-fsmooth" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
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
                            "This option needs an argument: "
                            "the sigma of the Gaussian "
                            "value to use." );
          nCurrentArg += 1;
        }

      } else if ( MATCH( sArg, "-revphaseflag" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-truncphaseflag" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-pial" )) {
        sPialName =  argv[nCurrentArg+1] ;
        nCurrentArg +=2 ;
      } else if ( MATCH( sArg, "-orig" )) {
        sOrigName =  argv[nCurrentArg+1] ;
        nCurrentArg +=2 ;
      } else if ( MATCH( sArg, "-segmentation-opacity" ) ||
		  MATCH( sArg, "-opacity" ) ||
                  MATCH( sArg, "-roialpha" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1 &&
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

      } else if ( MATCH( sArg, "-aux" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if (MATCH(sArg, "-aux-conform"  ) ) {

        /* set our flag */
        DebugNote( ("Enabling aux conform.") );
        bConformAuxVolume = TRUE;
        nCurrentArg ++;

      } else if ( MATCH( sArg, "-main-transform" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-aux-transform" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-label" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-headpts" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
             && '-' != argv[nCurrentArg+1][0] ) {

          /* read in the head pts and transform file name */
          DebugNote( ("Parsing -headpts option") );
          xUtil_strncpy( sHeadPts, argv[nCurrentArg+1], sizeof(sHeadPts) );

          /* if they gave us a transform file as well... */
          if ( argc > nCurrentArg + 2
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

        } else if ( MATCH( sArg, "-exit" ) ) {
          bExit = TRUE;
          nCurrentArg += 1;

        } else {

          /* misuse of that switch */
          tkm_DisplayError( "Parsing -headpts option",
                            "Expected one or two arguments",
                            "This option needs an argument: the file name "
                            "of the head points file, and optionally the "
                            "file name of the transform to use." );
          nCurrentArg += 1;
        }

      } else if ( MATCH( sArg, "-overlaycache" ) ) {

        /* check for the value following the switch */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-interface" ) ) {

        /* check for another value */
        if ( argc > nCurrentArg + 1
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

      } else if ( MATCH( sArg, "-" ) ) {

        /* set no edit mode. */
        DebugNote( ("Parsing - option") );
        bNoEdit = TRUE;
        nCurrentArg += 1;

      } else if ( MATCH( sArg, "-csurf" ) ) {

        /* set name of interface to tkmedit_csurf.tcl */
        DebugNote( ("Parsing -csurf option") );
        gbUseCsurfInterface = TRUE;
        nCurrentArg ++;

      } else if ( MATCH( sArg, "-mip" ) ) {

        /* Set the maximum intensity projection flag */
        DebugNote( ("Parsing -mip option") );
        bMIP = TRUE;
        nCurrentArg ++;

      } else {

        /* unrecognized option, build an error message and ignore it. */
        xUtil_snprintf( sError, sizeof(sError),
                        "Option %s not recognized", sArg );
        if (getenv("TK_EXIT_ON_CMD_ERROR")!=NULL) {
          printf( "\nParsing command line options %s.\n"
                  "This option was not recognized.\n", sError);
          exit(1);

        }
        tkm_DisplayError( "Parsing command line options",
                          sError,
                          "This option was not recognized and ignored." );
        nCurrentArg ++;
      }

      /* check for local keyword */
    } else if ( MATCH( sArg, "local" ) ||
                MATCH( sArg, "." ) ) {

      /* make sure subject is not already declared */
      if ( bSubjectDeclared ) {
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
         declared the subject yet. */
      if ( !bSubjectDeclared ) {

        /*make sure we have enough args and they arn't switches */
        if ( argc > nCurrentArg+1 ) {

          /* make sure the next two args arn't switches or local/. args */
          if ( '-' == argv[nCurrentArg+1][0] ) {

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

            if(!MATCH(sSubject,"getreg")){
              /* save subject home */
              DebugNote( ("Setting subject home from env") );
              eResult = SetSubjectHomeDirFromEnv( sSubject );
              DebugAssertThrow( (tkm_tErr_NoErr == eResult) );
	      //	    printf("Setting subject to %s\n",sSubject);
	      sprintf(tmpstr,"%s/%s",getenv("SUBJECTS_DIR"),sSubject);
	      if(!fio_FileExistsReadable(tmpstr)){
		printf("ERROR: cannot find subject %s in %s or it is not readable by you\n",
		       sSubject,getenv("SUBJECTS_DIR"));
		exit(1);
	      }
            } else {
              printf("Getting subject from registration file.\n");
              bGetSubjectFromReg = TRUE;
            }
            // Automatically set brightness/contrast for fsaverage
	    // Not needed as of 4/16/08
            if(0 && !strcmp(sSubject,"fsaverage")){
              fBrightnessMain = .58;
              fContrastMain = 14;
              bBrightContrastMain = TRUE;
            }
          }

          /* check for a surface. if we have enough args... */
          if ( argc > nCurrentArg ) {

            /* and this one isn't a switch or the local flag... */
            if ( '-' != argv[nCurrentArg][0]
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
  if ( bFatalError ) {
    /* probably have some error messages that didn't get printed cuz
       tcl wasn't loaded yet, so flush them to shell now. */
    PrintCachedTclErrorDlogsToShell();
    DebugPrint( ( "Fatal error in parsing command line args.\n" ) );
    exit(1);
  }

  /* if using local directory, copy 'local' in subject name for
     historical reasons. */
  if ( bLocalImageDir ) {
    DebugNote( ("Copying local into subject name") );
    xUtil_strncpy( sSubject, "local", sizeof(sSubject) );
    DebugNote( ("Copying empty string into image dir") );
    xUtil_strncpy( sImageDir, "", sizeof(sImageDir) );
  }

  /* disable editing */
  if ( bNoEdit ) {
    editflag = FALSE;
  }

  /* bleah. if we have the using mri read flag set, load the images from
     the subject name, otherwise use the image dir. */
  if ( bUsingMRIRead ) {
    DebugNote( ("Loading volume %s", sSubject) );
    eResult = LoadVolume( tkm_tVolumeType_Main, sSubject, bConformVolume );
    if ( tkm_tErr_NoErr != eResult ) {
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }
  } else {
    DebugNote( ("Loading volume %s", sImageDir) );
    eResult = LoadVolume( tkm_tVolumeType_Main, sImageDir, bConformVolume );
    if ( tkm_tErr_NoErr != eResult ) {
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }

    /* check to see if we don't have a subject */
    Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main],
                          sSubjectTest, sizeof(sSubjectTest) );
    if ( strcmp( sSubjectTest, "" ) == 0 ) {
      /* manually set the subject and image name */
      Volm_SetSubjectName( gAnatomicalVolume[tkm_tVolumeType_Main],
                           sSubject );
      Volm_SetVolumeName( gAnatomicalVolume[tkm_tVolumeType_Main],
                          sImageDir );
    }
  }

  /* If we're scaling up, do it now. Re-allocate the selection volume too. */
  if ( bScaleUpVolume ) {
    Volm_SetMinVoxelSizeToOne( gAnatomicalVolume[tkm_tVolumeType_Main] );
    Volm_GetMRIIdxToAnaIdxTransform( gAnatomicalVolume[tkm_tVolumeType_Main],
                                     &gMRIIdxToAnaIdxTransform );
    AllocateSelectionVolume();

    /* This changes when you resize the volume like that, so get it
       again. */
    Volm_GetMRIIdxToAnaIdxTransform( gAnatomicalVolume[tkm_tVolumeType_Main],
                                     &gMRIIdxToAnaIdxTransform );

  }

  /* If we got a non-default brightness and contrast or min and max,
     set it now. */
  if ( bBrightContrastMain ) {
    SetVolumeBrightnessAndContrast( tkm_tVolumeType_Main,
                                    fBrightnessMain, fContrastMain );
  }
  if ( bColorMain ) {
    SetVolumeColorMinMax( tkm_tVolumeType_Main, fColorMinMain, fColorMaxMain );
  }

  /* if reading in an aux image... */
  if ( bLoadingAuxVolume ) {
    DebugNote( ("Loading aux volume %s", sAuxVolume) );
    eResult = LoadVolume( tkm_tVolumeType_Aux, sAuxVolume, bConformAuxVolume );

    /* If we got a non-default brightness and contrast or min and max,
       set it now. */
    if ( bBrightContrastAux ) {
      SetVolumeBrightnessAndContrast( tkm_tVolumeType_Aux,
                                      fBrightnessAux, fContrastAux );
    }
    if ( bColorAux ) {
      SetVolumeColorMinMax( tkm_tVolumeType_Aux, fColorMinAux, fColorMaxAux );
    }
    if ( bScaleUpVolume ) {
      Volm_SetMinVoxelSizeToOne( gAnatomicalVolume[tkm_tVolumeType_Aux] );
      Volm_GetMRIIdxToAnaIdxTransform( gAnatomicalVolume[tkm_tVolumeType_Aux],
                                       &gMRIIdxToAnaIdxTransform );
    }
  }

  /* load in the display transforms. */
  if ( bLoadingMainTransform ) {
    DebugNote( ("Loading main display transform %s", sMainTransform) );
    eResult = LoadDisplayTransform( tkm_tVolumeType_Main, sMainTransform );
  }
  if ( bLoadingAuxTransform ) {
    DebugNote( ("Loading aux display transform %s", sAuxTransform) );
    eResult = LoadDisplayTransform( tkm_tVolumeType_Aux, sAuxTransform );
  }

  /* load surface. trasnsform must be inited first. */
  if ( bSurfaceDeclared ) {
    DebugNote( ("Loading surface") );
    eResult = LoadSurface( tkm_tSurfaceType_Main, sSurface );
  }

  /* load aux surface. */
  if ( bAuxSurfaceDeclared ) {
    DebugNote( ("Loading aux surface") );
    eResult = LoadSurface( tkm_tSurfaceType_Aux, sAuxSurface );
  }

  /* Import an annotation */
  if ( bLoadingAnnotation ) {

    if ( bSurfaceDeclared ) {
      DebugNote( ("Importing annotation") );
      ImportSurfaceAnnotationToSegmentation( tkm_tSegType_Main,
                                             sAnnotation,
                                             sAnnotationColorTable );
      /* set roi alpha */
      if ( bSegmentationAlpha ) {
        SetSegmentationAlpha( fSegmentationAlpha );
      }

    } else {
      tkm_DisplayError( "Importing Annotation",
                        "Surface not loaded",
                        "In order to import an annotation, you must also "
                        "load a surface." );
    }
  }

  /* load segmentation */
  if ( bLoadingSegmentation ) {
    eResult = LoadSegmentationVolume( tkm_tSegType_Main,
                                      sSegmentationPath,
                                      sSegmentationColorFile,
                                      bConformSegmentation );
    printf( "LoadSegmentationVolume main %s %s\n",
            sSegmentationPath, sSegmentationColorFile );
    /* set roi alpha */
    if ( bSegmentationAlpha ) {
      SetSegmentationAlpha( fSegmentationAlpha );
    }
  }

  /* load aux segmentation */
  if ( bLoadingAuxSegmentation ) {
    eResult = LoadSegmentationVolume( tkm_tSegType_Aux,
                                      sAuxSegmentationPath,
                                      sAuxSegmentationColorFile,
                                      bConformAuxSegmentation );
    printf( "LoadSegmentationVolume aux %s %s\n",
            sAuxSegmentationPath, sAuxSegmentationColorFile );
  }

  /* load VLIs */
  if ( bLoadingVLI ) {
    eResult = LoadVLIs( sVLIFile1, sVLIFile2 );
  }

  /* load the label */
  if ( bLoadingLabel ) {
    eResult = LoadSelectionFromLabelFile( sLabel );
  }

  /* load head pts */
  if ( bLoadingHeadPts ) {
    if ( bHaveHeadPtsTransform ) {
      eResult = LoadHeadPts( sHeadPts, sHeadPtsTransform );
    } else {
      eResult = LoadHeadPts( sHeadPts, NULL );
    }
  }

  /* load functional overlay data */
  if ( bLoadingOverlay ) {

    /* Check that they specified a reg type. */
    if ( FunD_tRegistration_None == overlayRegType ) {

      OutputPrint "INFO: No registration type specified for overlay, "
        "assuming identity.\n" EndOutputPrint;
      overlayRegType = FunD_tRegistration_Identity;
    }

    eResult = LoadFunctionalOverlay( sOverlayFileName,
                                     psOverlayOffsetFileName,
                                     overlayRegType,
                                     sOverlayRegistration );

    if ( eResult == tkm_tErr_NoErr && bSmooth ) {
      eResult = SmoothOverlayData( smoothSigma );
    }
  }

  /* load functional time course data */
  if ( bLoadingTimeCourse ) {

    /* Check that they specified a reg type. */
    if ( FunD_tRegistration_None == timecourseRegType ) {

      OutputPrint "INFO: No registration type specified for time course, "
        "assuming identity.\n" EndOutputPrint;
      timecourseRegType = FunD_tRegistration_Identity;
    }

    eResult = LoadFunctionalTimeCourse( sTimeCourseFileName,
                                        psTimeCourseOffsetFileName,
                                        timecourseRegType,
                                        sTimeCourseRegistration );
  }

  /* set registration */
  DebugNote( ("%sabling registration", bEnablingRegistration?"En":"Dis") );
  FunV_EnableRegistration( gFunctionalVolume, bEnablingRegistration );

  /* if regisration is enabled, set conversion method to round. */
  if ( bEnablingRegistration ) {
    FunV_SetConversionMethod( gFunctionalVolume,
                              FunD_tConversionMethod_Round );
  }

  if ( bThresh && bMax && !bSlope && !bMid ) {
    mid = min + (max-min)/2;
    bMid = 1;
    slope = 1/(max-min);
    bSlope = 1;
  }

  /* set functional color scale stuff */
  if ( bThresh || bMid || bSlope ) {
    DebugNote( ("Setting functional threshold min=%f mid=%f slope=%f\n",
                min, mid, slope ) );
    eFunctional = FunV_SetThreshold( gFunctionalVolume, min, mid, slope );
  }
  if ( bTruncPhaseFlag ) {
    DebugNote( ("Setting trunc phase flag to %d\n", nTruncPhaseFlag) );
    eFunctional = FunV_SetDisplayFlag( gFunctionalVolume,
                                       FunV_tDisplayFlag_Ol_TruncateNegative,
                                       (tBoolean) nTruncPhaseFlag );
  }
  if ( bRevPhaseFlag ) {
    DebugNote( ("Setting rev phase flag to %d\n", nRevPhaseFlag ) );
    eFunctional = FunV_SetDisplayFlag( gFunctionalVolume,
                                       FunV_tDisplayFlag_Ol_ReversePhase,
                                       (tBoolean) nRevPhaseFlag );
  }
  if ( bUseOverlayCacheFlag ) {
    DebugNote( ("Setting overlay cache flag to %d\n",
                nUseOverlayCacheFlag) );
    eFunctional = FunV_UseOverlayCache( gFunctionalVolume,
                                        (tBoolean) nUseOverlayCacheFlag );
  }

  /* Maximum intensity projection. */
  if ( bMIP ) {
    DebugNote( ("Setting MIP flag on\n") );
    MWin_SetDisplayFlag( gMeditWindow, -1,
                         DspA_tDisplayFlag_MaxIntProj, bMIP );
  }

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

  if ( NULL == gMeditWindow )
    return;

  MWin_GetWindowSize( gMeditWindow,
                      &xloc, &yloc, &width, &height );

  size = width*height;
  nNumBytes = sizeof( unsigned short ) * size;

  red  = (unsigned short*) malloc( nNumBytes );
  if ( NULL == red ) {
    DebugPrint( ( "save_rgb: allocation of red buffer failed.\n" ) );
    goto cleanup;
  }
  green = (unsigned short*) malloc( nNumBytes );
  if ( NULL == green ) {
    DebugPrint( ( "save_rgb: allocation of green buffer failed.\n" ) );
    goto cleanup;
  }
  blue  = (unsigned short*) malloc( nNumBytes );
  if ( NULL == blue ) {
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
  if ( GL_NO_ERROR != eGL ) {
    DebugPrint( ( "glReadPixels got error %d\n", eGL ) );
    goto cleanup;
  }
  glReadPixels(0, 0, width, height, GL_GREEN,
               GL_UNSIGNED_SHORT, (GLvoid *)green);
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL ) {
    DebugPrint( ( "glReadPixels got error %d\n", eGL ) );
    goto cleanup;
  }
  glReadPixels(0, 0, width, height, GL_BLUE,
               GL_UNSIGNED_SHORT, (GLvoid *)blue);
  eGL = glGetError ();
  if ( GL_NO_ERROR != eGL ) {
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
  if (fp==NULL) {
    printf("medit: ### can't create file %s\n",fname);
    return;
  }
  fclose(fp);
  image = iopen(fname,"w",UNCOMPRESSED(1), 3, width, height, 3);
  for (y = 0 ; y < height; y++) {
    r = red + y * width;
    g = green + y * width;
    b = blue + y * width;
    putrow(image, r, y, 0);
    putrow(image, g, y, 1);
    putrow(image, b, y, 2);
  }
  iclose(image);

 cleanup:

  if ( NULL != red )
    free( red );
  if ( NULL != green )
    free( green );
  if ( NULL != blue )
    free( blue );
}

void save_tiff (char* fname) {
  GLint rowlength, skiprows, skippixels, alignment;
  GLboolean swapbytes, lsbfirst;
  GLubyte* pixel_data = NULL;
  GLenum gl_error;
  TIFF *tiff = NULL;
  tsize_t line_bytes;
  unsigned char* line_buffer = NULL;
  int scan_line_size;
  int strip_size;
  int row;
  int x, y, height, width;

  if ( NULL == gMeditWindow )
    return;

  MWin_GetWindowSize( gMeditWindow,
                      &x, &y, &width, &height );

  /* Allocate a buffer for pixels. */
  pixel_data = (GLubyte*) malloc (width * height * 3);
  if (NULL == pixel_data) {
    DebugPrint(("Couldn't allocate pixel data."));
    goto error;
  }

  /* Read from the front buffer. */
  glReadBuffer (GL_FRONT);

  /* Save our unpack attributes. */
  glGetBooleanv (GL_PACK_SWAP_BYTES, &swapbytes);
  glGetBooleanv (GL_PACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv (GL_PACK_ROW_LENGTH, &rowlength);
  glGetIntegerv (GL_PACK_SKIP_ROWS, &skiprows);
  glGetIntegerv (GL_PACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv (GL_PACK_ALIGNMENT, &alignment);

  /* Set them. */
  glPixelStorei (GL_PACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei (GL_PACK_ROW_LENGTH, 0);
  glPixelStorei (GL_PACK_SKIP_ROWS, 0);
  glPixelStorei (GL_PACK_SKIP_PIXELS, 0);
  glPixelStorei (GL_PACK_ALIGNMENT, 1);

  /* Read RGB pixel data. */
  glReadPixels (0, 0, width, height, GL_RGB,
                GL_UNSIGNED_BYTE, (GLvoid*)pixel_data);

  /* Check error at this point. */
  gl_error = glGetError ();

  /* Restore the attributes. */
  glPixelStorei (GL_PACK_SWAP_BYTES, swapbytes);
  glPixelStorei (GL_PACK_LSB_FIRST, lsbfirst);
  glPixelStorei (GL_PACK_ROW_LENGTH, rowlength);
  glPixelStorei (GL_PACK_SKIP_ROWS, skiprows);
  glPixelStorei (GL_PACK_SKIP_PIXELS, skippixels);
  glPixelStorei (GL_PACK_ALIGNMENT, alignment);

  /* Handle error now. We can bail safely as we've restored the GL
     context. */
  if (GL_NO_ERROR != gl_error) {
    DebugPrint(("Error reading pixels."));
    goto error;
  }

  /* Open a TIFF. */
  tiff = TIFFOpen( fname, "w" );
  if (NULL == tiff) {
    DebugPrint(("Couldn't create file."));
    goto error;
  }

  /* Set the TIFF info. */
  TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

  /* Calculate some sizes and allocate a line buffer. */
  line_bytes = 3 * width;
  line_buffer = NULL;
  scan_line_size = TIFFScanlineSize (tiff);
  if (scan_line_size != line_bytes) {
    fprintf (stderr,"surfer: scan_line_size %d, line_bytes %d\n",
             scan_line_size, (int)line_bytes);
  }

  line_buffer = (unsigned char*) _TIFFmalloc( scan_line_size  );
  if (NULL == line_buffer) {
    DebugPrint(("Couldn't create tiff storage."));
    goto error;
  }

  /* Set the strip size to default. */
  strip_size = TIFFDefaultStripSize (tiff, width * 3);
  TIFFSetField (tiff, TIFFTAG_ROWSPERSTRIP, strip_size);

  /* Write line by line (bottom to top). */
  for (row = 0; row < height; row++) {
    memmove (line_buffer, &pixel_data[(height-row-1) * line_bytes],
             line_bytes);
    TIFFWriteScanline (tiff, line_buffer, row, 0);
  }

  goto cleanup;
 error:

 cleanup:

  if (NULL != pixel_data)
    free (pixel_data);

  /* Close the tiff file and free the line buffer. */
  if (NULL != tiff)
    TIFFClose (tiff);
  if (NULL != line_buffer)
    _TIFFfree (line_buffer);
}


void GotoSurfaceVertex ( Surf_tVertexSet iSurface, int inVertex ) {

  Surf_tErr eSurface = Surf_tErr_NoErr;
  MWin_tErr eWindow  = MWin_tErr_NoErr;
  xVoxel    anaIdx;
  char      sDescription[STRLEN];
  char      sSetName[STRLEN];

  /* make sure we have a surface. */
  if ( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the vertex */
  eSurface = Surf_GetNthVertex( gSurface[tkm_tSurfaceType_Main],
                                iSurface, inVertex, &anaIdx,
                                sDescription );
  if ( Surf_tErr_NoErr != eSurface )
    goto error;

  /* print the result string */
  Surf_GetSurfaceSetName( iSurface, sSetName );
  OutputPrint "%s vertex index %d:\n\t%s\n",
    sSetName, inVertex, sDescription EndOutputPrint;

  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if ( MWin_tErr_NoErr != eWindow )
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
  if ( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the cursor */
  eWindow = MWin_GetCursor ( gMeditWindow, &cursor );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  /* get the vertex */
  eSurface = Surf_GetClosestVertexVoxel( gSurface[tkm_tSurfaceType_Main],
                                         iSet, &cursor, &anaIdx,
                                         sDescription);
  if ( Surf_tErr_NoErr != eSurface )
    goto error;

  /* print the result string */
  Surf_GetSurfaceSetName( iSet, sSetName );
  OutputPrint "Nearest %s vertex to %d, %d, %d:\n\t%s\n",
    sSetName, xVoxl_ExpandInt( &cursor ), sDescription EndOutputPrint;

  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &anaIdx );
  if ( MWin_tErr_NoErr != eWindow )
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
  if ( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* get the cursor */
  eWindow = MWin_GetCursor ( gMeditWindow, &cursor );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  /* get the verteices */
  eWindow = MWin_GetClosestInterpSurfVoxel( gMeditWindow,
                                            tkm_tSurfaceType_Main,
                                            iSet, &cursor,
                                            &origAnaIdx, &interpAnaIdx,
                                            sDescription);
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  /* print the result string */
  Surf_GetSurfaceSetName( iSet, sSetName );
  OutputPrint "Nearest %s vertex to %d, %d, %d:\n\t%s\n",
    sSetName, xVoxl_ExpandInt( &cursor ), sDescription EndOutputPrint;

  /* tell the window to go there. */
  eWindow = MWin_SetCursor ( gMeditWindow, -1, &interpAnaIdx );
  if ( MWin_tErr_NoErr != eWindow )
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
  if ( NULL == gSurface[tkm_tSurfaceType_Main] )
    goto error;

  /* call the average function. */
  eSurface = Surf_AverageVertexPositions( gSurface[tkm_tSurfaceType_Main],
                                          inNumAverages );
  if ( Surf_tErr_NoErr != eSurface )
    goto error;

  /* redraw the window. */
  eWindow = MWin_RedrawAll( gMeditWindow );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  goto cleanup;

 error:

  DebugPrint( ( "Error in AverageSurfaceVertexPositions( %d )\n",
                (int)inNumAverages ) );

 cleanup:
  return;
}

void SetUseRealRAS ( tkm_tSurfaceType iType,
                     tBoolean ibUseRealRAS ) {


  gbUseRealRAS = ibUseRealRAS;
  tkm_SendTclCommand( tkm_tTclCommand_UpdateUseRealRAS,
                      gbUseRealRAS ? "1" : "0" );

  /* Recalc the surface transform. */
  CalcAndSetSurfaceClientTransformation( iType );

  /* Mark the surface dirty. */
  MWin_SetSurface( gMeditWindow, -1, iType, gSurface[iType] );

  /* redraw the window. */
  MWin_RedrawAll( gMeditWindow );
}



tkm_tErr CalcAndSetSurfaceClientTransformation ( tkm_tSurfaceType iType ) {

  tkm_tErr        eResult          = tkm_tErr_NoErr;
  Trns_tErr       eTrns            = Trns_tErr_NoErr;
  mriTransformRef surfaceTransform = NULL;
  MATRIX*         tmp1             = NULL;
  MATRIX*         tmp2             = NULL;

  DebugEnterFunction( ("CalcAndSetSurfaceClientTransformation( iType=%d )",
                       (int)iType) );

  // gMRIIdxToAnaIdxTransform keeps track of
  //        src ---> RAS
  //  non-   |        |
  //  triv   |        |identity
  //         V        V
  //    conformed --> RAS
  //      256^3
  DebugNote( ("Cloning gMRIIdxToAnaIdxTransform to get surfaceTransform") );
  eTrns = Trns_DeepClone( gMRIIdxToAnaIdxTransform, &surfaceTransform );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTrns),
                     eResult, tkm_tErr_CouldntAllocate );

  // surfaceTransform keeps track of
  //
  //           known
  //        src ---> RAS
  //         |        |
  // non-    |        | non-trivial
  // triv    V        V
  //   conformed---> RAS  (c_(r,a,s) = 0)
  //    256^3  standard
  //
  /* If we are using real RAS, we just use the extract_i_to_r
     transform for src->conformed, otherwise use the standard. */
  if ( gbUseRealRAS ) {

    tmp1 = MatrixCopy
      (extract_i_to_r
       ( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ),
       NULL );

  } else {

    Trns_GetBtoRAS(surfaceTransform, &tmp1);
    *MATRIX_RELT(tmp1, 1, 1) = -1;
    *MATRIX_RELT(tmp1, 2, 1) =  0;
    *MATRIX_RELT(tmp1, 3, 1) =  0;

    *MATRIX_RELT(tmp1, 1, 2) =  0;
    *MATRIX_RELT(tmp1, 2, 2) =  0;
    *MATRIX_RELT(tmp1, 3, 2) = -1;

    *MATRIX_RELT(tmp1, 1, 3) =  0;
    *MATRIX_RELT(tmp1, 2, 3) =  1;
    *MATRIX_RELT(tmp1, 3, 3) =  0;

    *MATRIX_RELT(tmp1, 1, 4) = 128;
    *MATRIX_RELT(tmp1, 2, 4) = -128;
    *MATRIX_RELT(tmp1, 3, 4) = 128;
  }

  Trns_CopyBtoRAS(surfaceTransform, tmp1);

  // in order to calculate ARAStoBRAS,
  // we use   ( RAS-> src ) then  ( src -> conformed )  then ( B->RAS )
  tmp1 = MatrixInverse(gMRIIdxToAnaIdxTransform->mAtoRAS, NULL);
  tmp2 = MatrixMultiply(gMRIIdxToAnaIdxTransform->mAtoB, tmp1, NULL);
  tmp1 = MatrixMultiply(surfaceTransform->mBtoRAS, tmp2, NULL);
  Trns_CopyARAStoBRAS(surfaceTransform, tmp1);

  Surf_SetTransform( gSurface[iType], surfaceTransform );

#if 0
  DebugPrint(("Set transform for surface %d\n", iType));
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

  goto cleanup;

 error:

 cleanup:

  if ( NULL != tmp1 )
    MatrixFree(&tmp1);
  if ( NULL != tmp2 )
    MatrixFree(&tmp2);

  DebugExitFunction;

  return eResult;
}


void CopyEditDatFileName ( char* osFileName, int izFileName ) {

  static tBoolean bWarnedLocal = FALSE;
  char*     pLocalEditFile  = NULL;
  tBoolean  bFoundLocal     = FALSE;
  char      sSubjectName[tkm_knNameLen] = "";
  FILE*     fTest           = NULL;
  char      sFileName[tkm_knPathLen] = "";

  /* First check if the local edit.dat file exists. If not, use the
     normal one. */
  bFoundLocal = FALSE;
  pLocalEditFile = getenv( "FS_SAVE_GOTO_POINT" );
  if ( NULL != pLocalEditFile ) {

    Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main],
                          sSubjectName, sizeof( sSubjectName ));
    sprintf( sFileName, "%s-%s", pLocalEditFile, sSubjectName );

    fTest = fopen( sFileName, "a" );
    if ( fTest ) {
      bFoundLocal = TRUE;
      fclose( fTest );

      if ( !bWarnedLocal ) {
        OutputPrint "tkmedit: Using local edit.dat file %s\n", sFileName
          EndOutputPrint;
        bWarnedLocal = TRUE;
      }
    }
  }

  if ( !bFoundLocal ) {

    /* Make the normal file name. */
    DebugNote( ("Making file name from edit.dat") );
    MakeFileName( "edit.dat", tkm_tFileName_Edit,
                  sFileName, sizeof(sFileName) );
  }

  /* Return the file name. */
  strncpy( osFileName, sFileName, izFileName );
}

void WriteVoxelToEditFile ( xVoxelRef iAnaIdx ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen] = "";
  FILE*     file            = NULL;
  xVoxel    MRIIdx;
  xVoxel    ras;
  xVoxel    tal;
  tBoolean  bHasTransform;

  DebugEnterFunction( ("WriteVoxelToEditFile ( iAnaIdx=%d,%d,%d )",
                       xVoxl_ExpandInt( iAnaIdx )) );

  /* Get the file name. */
  CopyEditDatFileName( sFileName, sizeof(sFileName) );

  /* open it */
  DebugNote( ("Opening edit file") );
  file = fopen( sFileName, "w" );
  DebugAssertThrowX( (NULL != file), eResult, tkm_tErr_ErrorAccessingFile );

  /* convert idx to ras */
  eVolume = Volm_ConvertIdxToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                     iAnaIdx, &MRIIdx );
  /* TESTING: ONLY WRITE SURFACE RAS */
  /*
    if ( gbUseRealRAS ) {
    eVolume = Volm_ConvertMRIIdxToRAS( gAnatomicalVolume[tkm_tVolumeType_Main],
    &MRIIdx, &ras );
    } else {
  */
  eVolume =
    Volm_ConvertMRIIdxToSurfaceRAS( gAnatomicalVolume[tkm_tVolumeType_Main],
                                    &MRIIdx, &ras );
  /*
    }
  */
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* write RAS space pt to file */
  DebugNote( ("Writing ras edit point %.2f,%.2f,%.2f to file",
              xVoxl_ExpandFloat( &ras ) ));
  fprintf( file,"%f %f %f\n", xVoxl_ExpandFloat( &ras ) );

  /* convert to tal and write that. */
  Volm_HasTalTransform( gAnatomicalVolume[tkm_tVolumeType_Main],
                        &bHasTransform );
  if ( bHasTransform ) {

    eVolume = Volm_ConvertIdxToTal( gAnatomicalVolume[tkm_tVolumeType_Main],
                                    iAnaIdx, &tal );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                       eResult, tkm_tErr_ErrorAccessingVolume );
    DebugNote( ("Writing tal edit point %.2f,%.2f,%.2f to file",
                xVoxl_ExpandFloat( &tal ) ));
    fprintf( file,"%f %f %f\n", xVoxl_ExpandFloat( &tal ) );
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  if ( NULL != file ) {
    DebugNote( ("Closing edit file") );
    fclose( file );
  }

  DebugExitFunction;
}

void GotoSavedCursor () {

  char sFileName[tkm_knPathLen] = "";

  DebugEnterFunction( ("GotoSavedCursor ()") );

  /* Get the file name. */
  CopyEditDatFileName( sFileName, sizeof(sFileName) );

  ReadCursorFromEditFile( sFileName );

  DebugExitFunction;
}

void GotoAnotherSubjectSavedCursor ( char* isSubjectName ) {}

void ReadCursorFromEditFile ( char* isFileName ) {

  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  FILE*      file    = NULL;
  float      fRASX   = 0;
  float      fRASY   = 0;
  float      fRASZ   = 0;
  xVoxel    ras;
  xVoxel    MRIIdx;
  xVoxel    idx;
  char      sError[tkm_knErrStringLen] = "";

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

  /* convert to screen voxel. */
  /* TESTING: ONLY WRITE SURFACE RAS */
  /*
    if ( gbUseRealRAS ) {
    eVolume = Volm_ConvertRASToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
    &ras, &MRIIdx );
    } else {
  */
  eVolume =
    Volm_ConvertSurfaceRASToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                    &ras, &MRIIdx );
  /*
    }
  */
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );
  eVolume = Volm_ConvertMRIIdxToIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                     &MRIIdx, &idx );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* build and set cursor */
  DebugNote( ("Setting cursor in main window") );
  MWin_SetCursor( gMeditWindow, -1, &idx );
  DebugNote( ("Setting zoom center in main window") );
  MWin_SetZoomCenterToCursor( gMeditWindow, -1 );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  sprintf( sError, "Tkmedit couldn't read the cursor from %s.", isFileName );
  tkm_DisplayError( "Going to saved cursor",
                    tkm_GetErrorString(eResult),
                    sError );
  EndDebugCatch;

  if ( NULL != file ) {
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

void tkm_HandleIdle () {}

// ================================================================== SURFACES

tkm_tErr LoadSurface ( tkm_tSurfaceType iType,
                       char*    isName ) {

  tkm_tErr  eResult           = tkm_tErr_NoErr;
  Surf_tErr eSurface          = Surf_tErr_NoErr;
  tBoolean  bLoaded           = FALSE;
  char      sName[tkm_knPathLen]       = "";
  char      sError[tkm_knErrStringLen] = "";
  tBoolean  bUseRealRAS;
  char      sTclArgs[tkm_knNameLen] = "";
  char*     sBaseName;
  char      sHemi[3];
  char      sFileNameCopy[tkm_knPathLen];

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
  if ( NULL != gSurface[iType] ) {
    Surf_Delete( &gSurface[iType] );
  }

  /* create the surface */
  DebugNote( ("Creating surface") );
  eSurface = Surf_New( &gSurface[iType], sName );
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
                     eResult, tkm_tErr_CouldntLoadSurface );

  /* see if it was loaded */
  DebugNote( ("Loading main vertex set") );
  eSurface = Surf_IsVertexSetLoaded( gSurface[iType], Surf_tVertexSet_Main,
                                     &bLoaded );
  DebugAssertThrowX( (bLoaded), eResult, tkm_tErr_CouldntLoadSurface );

  /* set the medit window surface. */
  DebugNote( ("Setting surface in main window") );
  MWin_SetSurface( gMeditWindow, -1, iType, gSurface[iType] );

  /* If the useRealRAS are the same or this is the first time we're
     setting it, call SetUseRealRAS automatically. If not, prompt the
     user which they want to use. */
  Surf_UsesRealRAS( gSurface[iType], &bUseRealRAS );
  if ( gbUseRealRAS == bUseRealRAS ||
       gbSetFirstUseRealRAS ) {
    SetUseRealRAS( iType, bUseRealRAS );
    gbSetFirstUseRealRAS = TRUE;
  } else if ( bUseRealRAS != gbUseRealRAS ) {
    tkm_SendTclCommand( tkm_tTclCommand_DoResolveUseRealRASDlog, "" );
  }

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

  /* Extract the hemisphere (lh or rh) and update in the tcl land. */
  DebugNote( ("Extracting hemi") );
  strncpy( sFileNameCopy, isName, sizeof(sFileNameCopy) );
  sBaseName = basename( sFileNameCopy );
  if ( strlen(sBaseName) > 2 &&
       (sBaseName[0] == 'l' || sBaseName[0] == 'r') &&
       sBaseName[1] == 'h' && sBaseName[2] == '.' ) {
    strncpy( sHemi, sBaseName, 2 );
    sprintf( sTclArgs, "%d %s", (int)iType, sHemi );
    tkm_SendTclCommand( tkm_tTclCommand_UpdateSurfaceHemi, sTclArgs );
  }

  /* load other vertex sets. See if they exist first. */
  if(LoadOrigSurf){
    DebugNote( ("Loading orig set") );
    LoadSurfaceVertexSet( iType, Surf_tVertexSet_Original, sOrigName );
  }
  if(LoadPialSurf){
    DebugNote( ("Loading pial set") );
    LoadSurfaceVertexSet( iType, Surf_tVertexSet_Pial, sPialName );
  }

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
  switch ( iSet ) {
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

  DebugEnterFunction( ("UnloadSurface( iType=%d )", (int)iType ) );

  DebugAssertThrowX( (NULL != gSurface[iType]),
                     eResult, tkm_tErr_InvalidParameter );
  if ( !gSurface[iType] )
    return;

  /* free the surface. */
  eSurface = Surf_Delete( &gSurface[iType] );
  if ( Surf_tErr_NoErr != eSurface ) {
    DebugPrint( ( "Surf error %d in UnloadSurface: %s\n",
                  eSurface, Surf_GetErrorString( eSurface ) ) );
  }

  /* if this is our main surface, disable our loading options */
  if ( tkm_tSurfaceType_Main == iType ) {
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

  DebugEnterFunction( ("WriteSurfaceValues( iType=%d, isFileName=%s )",
                       (int)iType, isFileName) );

  DebugAssertThrowX( (NULL != gSurface[iType]),
                     eResult, tkm_tErr_InvalidParameter );
  if ( !gSurface[iType] )
    return;

  /* Write the values for this surface. */
  Surf_WriteValues( gSurface[iType], isFileName );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

tkm_tErr LoadSurfaceAnnotation ( tkm_tSurfaceType iType,
                                 char*            isName ) {

  tkm_tErr      eResult  = tkm_tErr_NoErr;
  Surf_tErr     eSurface = Surf_tErr_NoErr;

  DebugEnterFunction( ("LoadSurfaceAnnotation( iType=%d, isName=%s )",
                       (int)iType, isName) );

  DebugAssertThrowX( (NULL != isName),
                     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (iType >= 0 && iType < tkm_knNumSurfaceTypes),
                     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != gSurface[iType]),
                     eResult, tkm_tErr_SurfaceNotLoaded );

  /* Load the annotation. */
  DebugNote( ("Loading annotation %s", isName) );
  eSurface = Surf_LoadAnnotation( gSurface[iType], isName );
  DebugAssertThrowX( (Surf_tErr_NoErr == eSurface),
                     eResult, tkm_tErr_CouldntLoadSurfaceAnnotation );

  /* Set the surface in the window to purge the cache */
  DebugNote( ("Setting surface in window") );
  MWin_SetSurface( gMeditWindow, -1, tkm_tSurfaceType_Main,
                   gSurface[tkm_tSurfaceType_Main] );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    ScaleOverlayRegisgtration ( atof(argv[1]), argv[2][0] );
  }

  return TCL_OK;
}

int TclLoadHeadPts ( ClientData inClientData, Tcl_Interp* inInterp,
                     int argc, char* argv[] ) {

  if ( argc != 3 ) {
    Tcl_SetResult
      ( inInterp,
        "wrong # args: LoadHeadPts points_file transform_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    WriteHeadPointsFile();
  }

  return TCL_OK;
}

int TclSetSelectedHeadPointLabel ( ClientData inClientData,
                                   Tcl_Interp* inInterp,
                                   int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetSelectedHeadPointLabel label",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    SetSelectedHeadPointLabel( argv[1] );
  }

  return TCL_OK;
}

int TclAlignSelectedHeadPointToMRIIdx ( ClientData inClientData,
                                        Tcl_Interp* inInterp,
                                        int argc, char* argv[] ) {

  xVoxel MRIIdx;

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: AlignSelectedHeadPointToMRIIdx",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    xVoxl_Set( &MRIIdx, atoi( argv[1] ), atoi( argv[2] ), atoi( argv[3] ) );
    AlignSelectedHeadPointToMRIIdx( &MRIIdx );
  }

  return TCL_OK;
}

int TclLoadGCA ( ClientData inClientData, Tcl_Interp* inInterp,
                 int argc, char* argv[] ) {

  if ( argc < 3 || argc > 4 ) {
    Tcl_SetResult
      ( inInterp,
        "wrong # args: LoadGCA volume_dir transform_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    LoadGCA ( argv[1], argv[2] );
    if (argc == 4) {
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
    Tcl_SetResult
      ( inInterp,
        "wrong # args: LoadGCARenormalization renorm_file",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    SaveGCA ( argv[1] );
  }

  return TCL_OK;
}

int TclReadVoxelLabels ( ClientData inClientData,
                         Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {
  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: ReadVoxelLabels vli1 vli2",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    save_rgb ( argv[1] );
  }

  return TCL_OK;
}

int TclSaveTIFF ( ClientData inClientData, Tcl_Interp* inInterp,
                  int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SaveTIFF filename:string",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    save_tiff ( argv[1] );
  }

  return TCL_OK;
}

int TclThresholdVolume ( ClientData inClientData, Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {

  if ( argc != 4 ) {
    Tcl_SetResult
      ( inInterp,
        "wrong # args: ThresholdVolume threshold_value "
        "{0=below,1=above} new_value",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    if ( atoi(argv[1]) ) {
      FlipVolume( mri_tOrientation_Sagittal );
    }
    if ( atoi(argv[2]) ) {
      FlipVolume( mri_tOrientation_Horizontal );
    }
    if ( atoi(argv[3]) ) {
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

  if ( gbAcceptingTclCommands ) {
    switch ( argv[1][0] ) {
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

  if ( gbAcceptingTclCommands ) {
    LoadVolume( tkm_tVolumeType_Main, argv[1], FALSE );
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

  if ( gbAcceptingTclCommands ) {
    LoadVolume( tkm_tVolumeType_Aux, argv[1], FALSE );

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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    SaveVolume( atoi( argv[1] ), NULL, FALSE );
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

  if ( gbAcceptingTclCommands ) {
    SaveVolume( atoi( argv[1] ), argv[2], TRUE );
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

  if ( gbAcceptingTclCommands ) {

    /* get a volume index */
    volume = atoi( argv[1] );

    /* make sure it's main or aux. if we have that volume, load the
       transform. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
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

  if ( gbAcceptingTclCommands ) {

    /* get a volume index */
    volume = atoi( argv[1] );

    /* make sure it's main or aux. if we have that volume, unload the
       transform. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
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
    Tcl_SetResult ( inInterp, "wrong # args: UnloadGCA", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    ClearUndoVolume();
  }

  return TCL_OK;
}

int TclSetVolumeColorScale ( ClientData inClientData, Tcl_Interp* inInterp,
                             int argc, char* argv[] ) {

  tkm_tVolumeType volume = tkm_tVolumeType_Main;

  if ( argc != 6 ) {
    Tcl_SetResult
      ( inInterp,
        "wrong # args: SetVolumeColorScale volume threshold:"
        "float squash:float min:float max:float",
        TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* get a volume index */
    volume = atoi( argv[1] );

    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
        SetVolumeBrightnessAndContrast( volume,
                                        atof( argv[2] ), atof( argv[3] ));
        SetVolumeColorMinMax( volume, atof( argv[4] ), atof( argv[5] ));
      }
    }
  }

  return TCL_OK;
}

int TclSetVolumeBrightnessContrast ( ClientData inClientData,
                                     Tcl_Interp* inInterp,
                                     int argc, char* argv[] ) {

  tkm_tVolumeType volume = tkm_tVolumeType_Main;

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp,
                    "wrong # args: SetVolumeBrightnessAndContrast volume "
                    "threshold:float squash:float", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* get a volume index */
    volume = atoi( argv[1] );

    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
        SetVolumeBrightnessAndContrast( volume,
                                        atof( argv[2] ), atof( argv[3] ));
      }
    }
  }

  return TCL_OK;
}

int TclSetVolumeMinMax ( ClientData inClientData, Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {

  tkm_tVolumeType volume = tkm_tVolumeType_Main;

  if ( argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeMinMax volume "
                    "min:float max:float", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* get a volume index */
    volume = atoi( argv[1] );

    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
        SetVolumeColorMinMax( volume, atof( argv[2] ), atof( argv[3] ));
      }
    }
  }

  return TCL_OK;
}

int TclSetVolumeSampleType ( ClientData inClientData, Tcl_Interp* inInterp,
                             int argc, char* argv[] ) {

  tkm_tVolumeType  volume = tkm_tVolumeType_Main;
  Volm_tSampleType type   = Volm_tSampleType_Nearest;

  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeSampleType volume "
                    "sampleType:0=nearest,1=trilinear,2=sinc",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* Get a volume index and sample type. */
    volume = atoi( argv[1] );
    type   = (Volm_tSampleType)atoi( argv[2] );

    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
        SetVolumeSampleType( volume, type );
      }
    }
  }

  return TCL_OK;
}

int TclSetVolumeResampleMethod ( ClientData inClientData, Tcl_Interp* inInterp,
                                 int argc, char* argv[] ) {

  tkm_tVolumeType      volume  = tkm_tVolumeType_Main;
  Volm_tResampleMethod method  = Volm_tResampleMethod_RAS;

  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetVolumeResampleMethod volume "
                    "method:0=RAS,1=slice", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* Get a volume index and sample type. */
    volume = atoi( argv[1] );
    method = (Volm_tResampleMethod)atoi( argv[2] );

    /* make sure it's main or aux. if we have that volume, set the brightness
       and contrast for it. */
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      if ( NULL != gAnatomicalVolume[ volume ] ) {
        SetVolumeResampleMethod( volume, method );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    compareType = atoi( argv[1] );
    if ( compareType > FunV_tFindStatsComp_Invalid &&
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

  if ( gbAcceptingTclCommands ) {
    ClearSelection ();
  }

  return TCL_OK;
}


int TclSetCursorToCenterOfVolume ( ClientData inClientData,
                                   Tcl_Interp* inInterp,
                                   int argc, char* argv[] ) {

  int volume;

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SetCursorToCenterOfVolume volume",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    volume = atoi( argv[1] );
    if ( volume == tkm_tVolumeType_Main || volume == tkm_tVolumeType_Aux ) {
      SetCursorToCenterOfVolume( volume );
    }
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    RestoreUndoList ();
  }

  return TCL_OK;
}

int TclNewControlPoint ( ClientData inClientData, Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {

  xVoxel cursor;
  xVoxel MRIIdx;

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: NewControlPoint",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    MWin_GetCursor ( gMeditWindow, &cursor );
    Volm_ConvertIdxToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                             &cursor, &MRIIdx );
    AddMRIIdxControlPoint( &MRIIdx, TRUE );
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

  if ( gbAcceptingTclCommands ) {
    WriteControlPointFile ();
  }

  return TCL_OK;
}

int TclLoadMainSurface ( ClientData inClientData, Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {

  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;

  switch ( argc ) {
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
                    "wrong # args: LoadMainSurface 0=main,1=aux "
                    "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    LoadSurface ( type, sFileName );
  }

  return TCL_OK;
}

int TclLoadPialSurface ( ClientData inClientData, Tcl_Interp* inInterp,
                         int argc, char* argv[] ) {

  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;

  switch ( argc ) {
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
                    "wrong # args: LoadPialSurface 0=main,1=aux "
                    "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( type, Surf_tVertexSet_Pial, sFileName );
  }

  return TCL_OK;
}

int TclLoadOriginalSurface ( ClientData inClientData, Tcl_Interp* inInterp,
                             int argc, char* argv[] ) {

  tkm_tSurfaceType type        = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;

  switch ( argc ) {
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
                    "wrong # args: LoadOriginalSurface 0=main,1=aux "
                    "surface_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    LoadSurfaceVertexSet( type, Surf_tVertexSet_Original, sFileName );
  }

  return TCL_OK;
}

int TclUnloadSurface ( ClientData inClientData, Tcl_Interp* inInterp,
                       int argc, char* argv[] ) {

  tkm_tSurfaceType type = tkm_tSurfaceType_Main;

  switch ( argc ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    UnloadSurface ( tkm_tSurfaceType_Main );
    UnloadSurface ( tkm_tSurfaceType_Aux );
  }

  return TCL_OK;
}

int TclLoadSurfaceAnnotation ( ClientData inClientData, Tcl_Interp* inInterp,
                               int argc, char* argv[] ) {

  tkm_tSurfaceType type = tkm_tSurfaceType_Main;

  if ( argc != 3 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadSurfaceAnnotation "
                    "0=main,1=aux file_name:string", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    type = atoi( argv[1] );
    LoadSurfaceAnnotation( type, argv[2] );
  }

  return TCL_OK;
}


int TclWriteSurfaceValues ( ClientData inClientData, Tcl_Interp* inInterp,
                            int argc, char* argv[] ) {

  tkm_tSurfaceType type       = tkm_tSurfaceType_Main;
  char*       sFileName = NULL;

  switch ( argc ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    AverageSurfaceVertexPositions ( atoi(argv[1]) );
  }

  return TCL_OK;
}


int TclSetUseRealRAS ( ClientData inClientData,
                       Tcl_Interp* inInterp,
                       int argc, char* argv[] ) {

  if ( argc != 2) {
    Tcl_SetResult ( inInterp, "wrong # args: SetUseRealRAS use_real_ras",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    SetUseRealRAS( tkm_tSurfaceType_Main, atoi(argv[1]) );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    LoadSegmentationVolume ( atoi(argv[1]), argv[2], argv[3], FALSE );
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

  if ( gbAcceptingTclCommands ) {
    RecomputeSegmentation( atoi(argv[1]) );
  }

  return TCL_OK;
}

int TclSelectSegLabelAtCursor ( ClientData inClientData,
                                Tcl_Interp* inInterp,
                                int argc, char* argv[] ) {

  if ( argc != 2 ) {
    Tcl_SetResult ( inInterp, "wrong # args: SelectSegLabelAtCursor "
                    "volume:0=main,1=aux", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    SelectSegLabelAtCursor( atoi(argv[1]) );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands &&
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

  if ( gbAcceptingTclCommands ) {

    /* if they gave us a filename, use it, else pass null. */
    if ( argc == 2 ) {
      SaveSegmentationVolume( atoi(argv[1]), NULL, FALSE );
    } else {
      SaveSegmentationVolume( atoi(argv[1]), argv[2], TRUE );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    SetSegmentationAlpha ( atof(argv[1]) );
  }

  return TCL_OK;
}

int TclLoadFunctionalOverlay ( ClientData inClientData,
                               Tcl_Interp* inInterp,
                               int argc, char* argv[] ) {

  char sFileName[tkm_knPathLen] = "";
  char* sRegistration = NULL;
  FunD_tRegistrationType regType = FunD_tRegistration_None;

  if ( argc != 2 && argc != 3 && argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadFunctionalOverlay "
                    "filename [registrationType [registration]]",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    xUtil_strncpy( sFileName, argv[1], sizeof(sFileName) );
    regType = FunD_tRegistration_Identity;

    if ( argc > 2 && argv[2][0] != '\0' ) {
      regType = (FunD_tRegistrationType) atoi( argv[2] );
    }

    if ( argc > 3 && argv[3][0] != '\0' ) {
      sRegistration = (char*) calloc( tkm_knPathLen, sizeof(char) );
      xUtil_strncpy( sRegistration, argv[3], tkm_knPathLen );
    }

    LoadFunctionalOverlay( sFileName, NULL, regType, sRegistration );

    if ( NULL != sRegistration ) {
      free( sRegistration );
    }
  }

  return TCL_OK;
}

int TclLoadFunctionalTimeCourse ( ClientData inClientData,
                                  Tcl_Interp* inInterp,
                                  int argc, char* argv[] ) {

  char sFileName[tkm_knPathLen] = "";
  char* sRegistration = NULL;
  FunD_tRegistrationType regType = FunD_tRegistration_None;

  if ( argc != 2 && argc != 3 && argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadFunctionalTimeCourse "
                    "filename [registrationType [registration]]",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    xUtil_strncpy( sFileName, argv[1], sizeof(sFileName) );
    regType = FunD_tRegistration_Identity;

    if ( argc > 2 && argv[2][0] != '\0' ) {
      regType = (FunD_tRegistrationType) atoi( argv[2] );
    }

    if ( argc > 3 && argv[3][0] != '\0' ) {
      sRegistration = (char*) calloc( tkm_knPathLen, sizeof(char) );
      xUtil_strncpy( sRegistration, argv[3], tkm_knPathLen );
    }

    LoadFunctionalTimeCourse( sFileName, NULL, regType, sRegistration );

    if ( NULL != sRegistration ) {
      free( sRegistration );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    SetDTIAlpha ( atof(argv[1]) );
  }

  return TCL_OK;
}

int TclLoadGDF ( ClientData inClientData, Tcl_Interp* inInterp,
                 int argc, char* argv[] ) {

  tkm_tErr eResult = tkm_tErr_NoErr;

  if ( argc != 3 && argc != 4 ) {
    Tcl_SetResult ( inInterp, "wrong # args: LoadGDF filename:string"
                    "registrationType:int [registration:string]",
                    TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {
    if ( 3 == argc ) {
      eResult = LoadGDFHeader( argv[1], atoi(argv[2]), NULL );
    } else {
      eResult = LoadGDFHeader( argv[1], atoi(argv[2]), argv[3] );
    }
  }

  return TCL_OK;
}

int TclGetGDFID ( ClientData inClientData, Tcl_Interp* inInterp,
                  int argc, char* argv[] ) {

  if ( argc != 1 ) {
    Tcl_SetResult ( inInterp, "wrong # args: GetGDFID", TCL_VOLATILE );
    return TCL_ERROR;
  }

  if ( gbAcceptingTclCommands ) {

    /* Return the GDF ID. */
    Tcl_SetObjResult( inInterp, Tcl_NewIntObj( gGDFID ) );
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    bStatus = (tBoolean) atoi( argv[1] );
    if ( bStatus && !gbTimerOn ) {
      StartTimer();
    } else if ( !bStatus && gbTimerOn ) {
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

  if ( gbAcceptingTclCommands ) {
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

  if ( gbAcceptingTclCommands ) {
    if ( gMeditWindow != NULL ) {
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

  if ( gbAcceptingTclCommands ) {

    /* get the volume type and check it */
    volume = (tkm_tVolumeType) atoi(argv[1]);
    if ( volume != tkm_tVolumeType_Main && volume != tkm_tVolumeType_Aux ) {
      Tcl_SetResult ( inInterp, "wrong # args: GetSubjectName 0=main|1=aux",
                      TCL_VOLATILE );
      return TCL_ERROR;
    }

    /* get the subject name */
    eVolume = Volm_CopySubjectName( gAnatomicalVolume[volume], sSubjectName,
                                    sizeof( sSubjectName ));
    if ( Volm_tErr_NoErr != eVolume ) {
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

  if ( gbAcceptingTclCommands ) {

    /* get the volume type and check it */
    volume = (tkm_tVolumeType) atoi(argv[1]);
    if ( volume != tkm_tVolumeType_Main && volume != tkm_tVolumeType_Aux ) {
      Tcl_SetResult ( inInterp, "wrong # args: GetSubjectDir 0=main|1=aux",
                      TCL_VOLATILE );
      return TCL_ERROR;
    }

    /* get the source dir */
    eVolume = Volm_CopySourceDir( gAnatomicalVolume[volume], sSubjectDir,
                                  sizeof( sSubjectDir ));
    if ( Volm_tErr_NoErr != eVolume ) {
      Tcl_SetResult
        ( inInterp,
          "Couldn't get volume directory. Is volume loaded?",
          TCL_VOLATILE );
      return TCL_ERROR;
    }

    /* return it */
    Tcl_SetResult( inInterp, sSubjectDir, TCL_VOLATILE );
  }

  return TCL_OK;
}


/*=======================================================================*/

/* for tcl/tk */
#ifndef Windows_NT
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
#endif // Windows_NT
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
//static Tk_Window mainWindow;
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

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

  /* init our debugging macro code, if any. */
  InitDebugging( "tkmedit" );
  EnableDebuggingOutput;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  //printf("Gdiag_no = %d\n",Gdiag_no);

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
  for ( nArg = 0; nArg < argc; nArg++ ) {
    DebugPrint( ( "%s ", argv[nArg] ) );
  }
  DebugPrint( ( "\n\n" ) );
  DebugPrint
    (
      (
        "$Id: tkmedit.c,v 1.346 2011/10/15 15:43:39 fischl Exp $ $Name:  $\n"
        )
      );

  /* init glut */
  DebugNote( ("Initializing glut") );
  glutInit( &argc, argv );
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );

  /* init main window */
  DebugNote( ("Creating main window") );
  eWindow = MWin_New( &gMeditWindow, "", 512, 512 );
  DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );

  /* init our control pt list */
  DebugNote( ("Initializing control point list") );
  eResult = InitControlPointList();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );

  /* init other things */
  DebugNote( ("Initalizing undo list") );
  eResult = InitUndoList();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );

  DebugNote( ("Initalizing selection module") );
  eResult = InitSelectionModule();
  DebugAssertThrow( (eResult == tkm_tErr_NoErr) );

  DebugNote( ("Initalizing user cancel listener") );
  xUtil_InitializeUserCancel();

  /* create functional volume */
  DebugNote( ("Creating functional volume") );
  eFunctional = FunV_New( &gFunctionalVolume,
                          UpdateAndRedraw,
                          tkm_SendTclCommand,
                          SendTCLCommand );
  DebugAssertThrow( (FunV_tErr_NoError == eFunctional) );

  /* set windows data sources */
  DebugNote( ("Setting control points space in window") );
  eWindow = MWin_SetControlPointsSpace( gMeditWindow, -1, gControlPointList );
  DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );

  DebugNote( ("Setting selection list in window.") );
  eWindow = MWin_SetSelectionSpace( gMeditWindow, -1, gSelectionVolume );
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
  tkm_SendTclCommand( tkm_tTclCommand_ShowGDFOptions, "0" );

  /* init variables */
  for ( volume = 0; volume < tkm_knNumVolumeTypes; volume++ ) {
    gAnatomicalVolume[volume] = NULL;
    SetVolumeDirty( volume, FALSE );
  }

  for ( volume = 0; volume < tkm_knNumSegTypes; volume++ ) {
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
     in FREESURFER_HOME/tktools. if we can't find it, exit. */
  pFile = NULL;
  bFoundInterface = FALSE;

  /* did they request one? */
  if ( !MATCH( gInterfaceScriptName, "" ) ) {

    /* try to open it. */
    xUtil_strncpy( sInterfaceFileName, gInterfaceScriptName,
                   sizeof(sInterfaceFileName) );
    DebugNote( ( "Trying to open %s\n", sInterfaceFileName ) );
    pFile = fopen ( sInterfaceFileName,"r" );
    if ( NULL != pFile ) {
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

  if ( !bFoundInterface ) {

    /* next, look in TKMEDIT_SCRIPTS_DIR */
    DebugNote( ("Getting TKMEDIT_SCRIPTS_DIR env var") );
    pEnvVar = getenv("TKMEDIT_SCRIPTS_DIR");
    if ( NULL != pEnvVar) {
      xUtil_snprintf( sInterfaceFileName, sizeof(sInterfaceFileName),
                      "%s/tkmedit.tcl", pEnvVar);
      DebugNote( ( "Trying to open %s", sInterfaceFileName ) );
      pFile = fopen( sInterfaceFileName,"r" );
      if ( NULL != pFile ) {
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


  if ( !bFoundInterface ) {

    /* look in the local directory */
    DebugNote( ("Trying local tkmedit.tcl") );
    xUtil_strncpy( sInterfaceFileName, "tkmedit.tcl",
                   sizeof(sInterfaceFileName) );
    pFile = fopen( sInterfaceFileName, "r" );
    if ( NULL != pFile ) {
      bFoundInterface = TRUE;
      fclose( pFile );
    }
  }

  if ( !bFoundInterface ) {

    /* finally try in FREESURFER_HOME/tktools.
       make sure we have FREESURFER_HOME
       defined */
    DebugNote( ("Getting FREESURFER_HOME env var") );
    pEnvVar = getenv("FREESURFER_HOME");
    if ( NULL == pEnvVar) {
      tkm_DisplayError
        ( "Trying to find interface file",
          "No valid file found",
          "Tkmedit couldn't find a valid interface file "
          "(tkmedit.tcl). Normally this is in the directory "
          "specified in the FREESURFER_HOME varible, but this was "
          "not set in your environment. Tkmedit needs this "
          "file to run." );
      PrintCachedTclErrorDlogsToShell();
      exit( 1 );
    }

    xUtil_snprintf( sInterfaceFileName, sizeof(sInterfaceFileName),
                    "%s/tktools/%s", pEnvVar, "tkmedit.tcl");
    DebugNote( ( "Trying to open %s", sInterfaceFileName ) );
    pFile = fopen( sInterfaceFileName,"r" );
    if ( NULL != pFile ) {
      bFoundInterface = TRUE;
      fclose( pFile );
    }
  }

  /* if file still not found bail out. */
  if ( !bFoundInterface ) {
    tkm_DisplayError
      ( "Trying to find interface file",
        "No valid file found",
        "Tkmedit couldn't find a valid interface file. "
        "Normally this is in the directory specified in "
        "the FREESURFER_HOME varible, but it can also be in the same "
        "directory as tkmedit. Tkmedit needs this file to "
        "run. Try reinstalling your distribution." );
    PrintCachedTclErrorDlogsToShell();
    exit( 1 );
  }

  DebugPrint( ( "Using interface file %s\n", sInterfaceFileName) );

  /* process ctrl pts file */
  DebugNote( ("Reading control point file.") );
  ReadControlPointFile();

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

  /* start tcl/tk; first make interpreter */
  DebugNote( ("Creating Tcl interpreter") );
  interp = Tcl_CreateInterp();

  // kt - set global interp
  DebugNote( ("Setting global interpreter") );
  SetTclInterp ( interp );

  /* read tcl/tk internal startup scripts */
  eTcl = Tcl_Init( interp );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Tcl_Init returned %d: %s\n", (int)eTcl, interp->result) );
    tkm_DisplayError( "Initializing Tcl",
                      "Error initializing Tcl",
                      "For some reason, Tcl couldn't be initialized. Possible "
                      "reasons include it not being installed, installed "
                      "incorrectly, or the TCL_LIBRARY environment variable "
                      "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTcl );
  }
  eTcl = Tk_Init( interp );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Tk_Init returned %d: %s\n", (int)eTcl, interp->result) );
    tkm_DisplayError( "Initializing Tk",
                      "Error initializing Tk",
                      "For some reason, Tk couldn't be initialized. Possible "
                      "reasons include it not being installed, installed "
                      "incorrectly, or the TK_LIBRARY environment variable "
                      "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTk );
  }

  ////////// do the following only for RedHat Enterprise Linux
#if NEEDS_ITCL_ITK
  eTcl = Itcl_Init(interp);
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Itlc_Init returned %d: %s\n", (int)eTcl, interp->result) );
  }
  eTcl = Itk_Init(interp);
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Itk_Init returned %d: %s\n", (int)eTcl, interp->result) );
  }
#endif
  ///////////////////////////////////////////////////////////////////
  eTcl = Tix_Init( interp );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Tix_Init returned %d: %s\n", (int)eTcl, interp->result) );
    tkm_DisplayError( "Initializing Tix",
                      "Error initializing Tix",
                      "For some reason, Tix couldn't be initialized. Possible "
                      "reasons include it not being installed, installed "
                      "incorrectly, or the TIX_LIBRARY environment variable "
                      "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitTix );
  }
  eTcl = Blt_Init( interp );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Blt_Init returned %d: %s\n", (int)eTcl, interp->result) );
    tkm_DisplayError( "Initializing BLT",
                      "Error initializing BLT",
                      "For some reason, BLT couldn't be initialized. Possible "
                      "reasons include it not being installed, installed "
                      "incorrectly, or the TIX_LIBRARY environment variable "
                      "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitBLT );
  }

  Tcl_StaticPackage( interp, "BLT", Blt_Init, Blt_SafeInit );
#ifdef Windows_NT
#define Tix_SafeInit NULL
#endif
  Tcl_StaticPackage( interp, "Tix", Tix_Init, Tix_SafeInit );

  /* Initialize our Fsgdf functions. This is in fsgdf_wrap.c */
  eTcl = Fsgdf_Init( interp );
  if ( TCL_OK != eTcl ) {
    DebugPrint( ("Fsgdf_Init returned %d: %s\n", (int)eTcl, interp->result) );
    tkm_DisplayError( "Initializing FSGDF",
                      "Error initializing FGSDF",
                      "For some reason, the FSGDF code couldn't be "
                      "initialized. Possible reasons include it not being "
                      "installed, installed incorrectly, or the "
                      "FREESURFER_HOME environment variable "
                      "not being set or being set incorrectly." );
    DebugAssertThrowX( (TCL_OK == eTcl), eResult, tkm_tErr_CouldntInitBLT );
  }


#if SET_TCL_ENV_VAR
  /* restore env vars */
  if ( bChangedEnvVar ) {
    sprintf( sTclEnvVar, "%s=%s", "TCL_LIBRARY", sSavedTclLib );
    if ( putenv( sTclEnvVar ) ) {
      OutputPrint "ERROR: Couldn't restore TCL_LIBRARY env var.\n"
        EndOutputPrint;
      exit( 1 );
    }
    sprintf( sTkEnvVar, "%s=%s", "TK_LIBRARY", sSavedTkLib );
    if ( putenv( sTkEnvVar ) ) {
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
  tty    = isatty(0);
  DebugNote( ("Setting tcl_interactive var to %d", (int)tty) );
  Tcl_SetVar( interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY );

  /* register tcl commands */
  DebugNote( ("Registering tkmedit tcl commands") );

  Tcl_CreateCommand ( interp, "SetGCADisplayStatus",
                      (Tcl_CmdProc*) TclSetGCADisplayStatus,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "RecomputeSegmentation",
                      (Tcl_CmdProc*) TclRecomputeSegmentation,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "SelectSegLabelAtCursor",
                      (Tcl_CmdProc*) TclSelectSegLabelAtCursor,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "RestorePreviousSegmentation",
                      (Tcl_CmdProc*) TclRestorePreviousSegmentation,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "CrashHard",
                      (Tcl_CmdProc*) TclCrashHard,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ReadVoxelLabels",
                      (Tcl_CmdProc*) TclReadVoxelLabels,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "DebugPrint",
                      (Tcl_CmdProc*) TclDebugPrint,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RotateOverlayRegistration",
                      (Tcl_CmdProc*) TclRotateOverlayRegistration,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "TranslateOverlayRegistration",
                      (Tcl_CmdProc*) TclTranslateOverlayRegistration,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ScaleOverlayRegistration",
                      (Tcl_CmdProc*) TclScaleOverlayRegistration,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadHeadPts",
                      (Tcl_CmdProc*) TclLoadHeadPts,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RotateHeadPts",
                      (Tcl_CmdProc*) TclRotateHeadPts,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "TranslateHeadPts",
                      (Tcl_CmdProc*) TclTranslateHeadPts,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RestoreHeadPts",
                      (Tcl_CmdProc*) TclRestoreHeadPts,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "WriteHeadPointsTransform",
                      (Tcl_CmdProc*) TclWriteHeadPointsTransform,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "WriteHeadPointsFile",
                      (Tcl_CmdProc*) TclWriteHeadPointsFile,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetSelectedHeadPointLabel",
                      (Tcl_CmdProc*) TclSetSelectedHeadPointLabel,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "AlignSelectedHeadPointToMRIIdx",
                      (Tcl_CmdProc*) TclAlignSelectedHeadPointToMRIIdx,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadGCA",
                      (Tcl_CmdProc*) TclLoadGCA,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadGCARenorm",
                      (Tcl_CmdProc*) TclLoadGCARenormalization,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveGCA",
                      (Tcl_CmdProc*) TclSaveGCA,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveRGB",
                      (Tcl_CmdProc*) TclSaveRGB,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveTIFF",
                      (Tcl_CmdProc*) TclSaveTIFF,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ThresholdVolume",
                      (Tcl_CmdProc*) TclThresholdVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "FlipVolume",
                      (Tcl_CmdProc*) TclFlipVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RotateVolume",
                      (Tcl_CmdProc*) TclRotateVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadVolume",
                      (Tcl_CmdProc*) TclLoadVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadAuxVolume",
                      (Tcl_CmdProc*) TclLoadAuxVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UnloadVolume",
                      (Tcl_CmdProc*) TclUnloadVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveVolume",
                      (Tcl_CmdProc*) TclSaveVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveVolumeAs",
                      (Tcl_CmdProc*) TclSaveVolumeAs,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadVolumeDisplayTransform",
                      (Tcl_CmdProc*) TclLoadVolumeDisplayTransform,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UnloadVolumeDisplayTransform",
                      (Tcl_CmdProc*) TclUnloadVolumeDisplayTransform,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UnloadGCA",
                      (Tcl_CmdProc*) TclUnloadGCA,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SnapshotVolume",
                      (Tcl_CmdProc*) TclSnapshotVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RestoreVolumeFromSnapshot",
                      (Tcl_CmdProc*) TclRestoreVolumeFromSnapshot,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ClearUndoVolume",
                      (Tcl_CmdProc*) TclClearUndoVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetVolumeColorScale",
                      (Tcl_CmdProc*) TclSetVolumeColorScale,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetVolumeBrightnessContrast",
                      (Tcl_CmdProc*) TclSetVolumeBrightnessContrast,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetVolumeMinMax",
                      (Tcl_CmdProc*) TclSetVolumeMinMax,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetVolumeSampleType",
                      (Tcl_CmdProc*) TclSetVolumeSampleType,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetVolumeResampleMethod",
                      (Tcl_CmdProc*) TclSetVolumeResampleMethod,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveLabel",
                      (Tcl_CmdProc*) TclSaveLabel,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadLabel",
                      (Tcl_CmdProc*) TclLoadLabel,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "GraphSelectedRegion",
                      (Tcl_CmdProc*) TclGraphSelectedRegion,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SelectVoxelsByFuncValue",
                      (Tcl_CmdProc*) TclSelectVoxelsByFuncValue,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ClearSelection",
                      (Tcl_CmdProc*) TclClearSelection,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UndoLastEdit",
                      (Tcl_CmdProc*) TclUndoLastEdit,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetCursorToCenterOfVolume",
                      (Tcl_CmdProc*) TclSetCursorToCenterOfVolume,
                      (ClientData) 1, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SendCursor",
                      (Tcl_CmdProc*) TclSendCursor,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ReadCursor",
                      (Tcl_CmdProc*) TclReadCursor,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "NewControlPoint",
                      (Tcl_CmdProc*) TclNewControlPoint,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "WriteControlPointFile",
                      (Tcl_CmdProc*) TclWriteControlPointFile,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadMainSurface",
                      (Tcl_CmdProc*) TclLoadMainSurface,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadPialSurface",
                      (Tcl_CmdProc*) TclLoadPialSurface,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadOriginalSurface",
                      (Tcl_CmdProc*) TclLoadOriginalSurface,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UnloadSurface",
                      (Tcl_CmdProc*) TclUnloadSurface,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadSurfaceAnnotation",
                      (Tcl_CmdProc*) TclLoadSurfaceAnnotation,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "UnloadAllSurfaces",
                      (Tcl_CmdProc*) TclUnloadAllSurfaces,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "WriteSurfaceValues",
                      (Tcl_CmdProc*) TclWriteSurfaceValues,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "GotoMainVertex",
                      (Tcl_CmdProc*) TclGotoMainVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "GotoPialVertex",
                      (Tcl_CmdProc*) TclGotoPialVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "GotoOriginalVertex",
                      (Tcl_CmdProc*) TclGotoOriginalVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestMainVertex",
                      (Tcl_CmdProc*) TclShowNearestMainVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestOriginalVertex",
                      (Tcl_CmdProc*) TclShowNearestOriginalVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestPialVertex",
                      (Tcl_CmdProc*) TclShowNearestPialVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedMainVertex",
                      (Tcl_CmdProc*) TclShowNearestInterpolatedMainVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedOriginalVertex",
                      (Tcl_CmdProc*) TclShowNearestInterpolatedOriginalVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ShowNearestInterpolatedPialVertex",
                      (Tcl_CmdProc*) TclShowNearestInterpolatedPialVertex,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "AverageSurfaceVertexPositions",
                      (Tcl_CmdProc*) TclAverageSurfaceVertexPositions,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetUseRealRAS",
                      (Tcl_CmdProc*) TclSetUseRealRAS,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "NewSegmentationVolume",
                      (Tcl_CmdProc*) TclNewSegmentationVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand ( interp, "LoadSegmentationVolume",
                      (Tcl_CmdProc*) TclLoadSegmentationVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SaveSegmentationVolume",
                      (Tcl_CmdProc*) TclSaveSegmentationVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ExportChangedSegmentationVolume",
                      (Tcl_CmdProc*) TclExportChangedSegmentationVolume,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "ImportSurfaceAnnotationToSegmentation",
                      (Tcl_CmdProc*) TclImportSurfaceAnnotationToSegmentation,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetSegmentationAlpha",
                      (Tcl_CmdProc*) TclSetSegmentationAlpha,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadFunctionalOverlay",
                      (Tcl_CmdProc*) TclLoadFunctionalOverlay,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadFunctionalTimeCourse",
                      (Tcl_CmdProc*) TclLoadFunctionalTimeCourse,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadDTIVolumes",
                      (Tcl_CmdProc*) TclLoadDTIVolumes,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetDTIAlpha",
                      (Tcl_CmdProc*) TclSetDTIAlpha,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "LoadGDF",
                      (Tcl_CmdProc*) TclLoadGDF,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "GetGDFID",
                      (Tcl_CmdProc*) TclGetGDFID,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SmoothFunctionalOverlay",
                      (Tcl_CmdProc*) TclSmoothFunctionalOverlay,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "SetTimerStatus",
                      (Tcl_CmdProc*) TclSetTimerStatus,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "QuitMedit",
                      (Tcl_CmdProc*) TclQuitMedit,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand ( interp, "RedrawScreen",
                      (Tcl_CmdProc*) TclRedrawScreen,
                      (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand( interp, "StartTimer",
                     (Tcl_CmdProc*) TclStartTimer,
                     (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand( interp, "StopTimer",
                     (Tcl_CmdProc*) TclStopTimer,
                     (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand( interp, "ExecuteQueuedScripts",
                     (Tcl_CmdProc*) TclExecuteQueuedScripts,
                     (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  Tcl_CreateCommand( interp, "GetSubjectName",
                     (Tcl_CmdProc*) TclGetSubjectName,
                     (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );
  Tcl_CreateCommand( interp, "GetSubjectDir",
                     (Tcl_CmdProc*) TclGetSubjectDir,
                     (ClientData) 0, (Tcl_CmdDeleteProc*) NULL );

  /* parse the interface file */
  DebugNote( ("Parsing %s", sInterfaceFileName) );
  eTcl = Tcl_EvalFile( interp, sInterfaceFileName );
  if ( *interp->result != 0 ) {
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
#ifndef Windows_NT
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
#endif // Windows_NT
  if (tty) Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /* finish initializing functional volume */
  DebugNote( ("Initializing graph window") );
  eFunctional = FunV_InitGraphWindow( gFunctionalVolume, interp );
  if ( FunV_tErr_NoError != eFunctional ) {
    exit( 1 );
  }

  DebugNote( ("Registering functional Tcl commands") );
  eFunctional = FunV_RegisterTclCommands( gFunctionalVolume, interp );
  if ( FunV_tErr_NoError != eFunctional ) {
    exit( 1 );
  }

  /* set data in window */
  DebugNote( ("Setting functional volume in main window") );
  MWin_SetOverlayVolume( gMeditWindow, -1, gFunctionalVolume );

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
  SendVolumeColorScaleUpdate( tkm_tVolumeType_Main );
  if ( NULL != gAnatomicalVolume[tkm_tVolumeType_Aux] ) {
    SendVolumeColorScaleUpdate( tkm_tVolumeType_Aux );
  }

  /* set default segmentation alpha */
  DebugNote( ("Setting default segmentation alpha") );
  SetSegmentationAlpha( gfSegmentationAlpha );

  /* set default DTI alpha */
  DebugNote( ("Setting default DTI alpha") );
  SetDTIAlpha( gfDTIAlpha );

  /* if using csurf interface, call the func that hides a bunch of stuff. */
  if ( gbUseCsurfInterface ) {
    DebugNote( ("Enabling csurf version of interface") );
    tkm_SendTclCommand( tkm_tTclCommand_CsurfInterface, "" );
  }

  /* if the tktimer autostart environment variable is declared,
     automatically start the timer here. */
  if ( getenv( ksTkTimerAutoStartEnvVar ) ) {
    StartTimer();
  }

  /* never returns */
  if (bExit) {
    DebugNote( ("exiting as user requests") );
    exit(0) ;
  }
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
  if ( getenv( ksTkTimerAutoStartEnvVar ) ) {
    StopTimer();
  }
  /* delete window */
  DebugNote( ( "Deleting main window.\n" ) );
  MWin_Delete( &gMeditWindow );

  /* delete everything we allocated before */
  if ( NULL != gGCAVolume ) {
    DebugNote( ("Deleting GCA.") );
    DeleteGCA();
  }

  if ( NULL != gVLI1 || NULL != gVLI2 ) {
    DebugNote( ("Deleting VLIs.") );
    DeleteVLIs();
  }

  DebugNote( ("Deleteing functional volume.") );
  FunV_Delete( &gFunctionalVolume );

  DebugNote( ("Deleting selection module.") );
  DeleteSelectionModule ();

  if ( NULL != gControlPointList ) {
    DebugNote( ("Deleting control point list.") );
    x3Lst_Delete( &gControlPointList );
  }

  DebugNote( ("Deleting undo volume.") );
  DeleteUndoVolume ();

  DebugNote( ("Deleting undo list.") );
  DeleteUndoList ();

  DebugPrint( ("Program terminated successfully.\n") );
  DebugNote( ( "Deleting debugging and exiting.\n" ) );
  DeleteDebugging;

  /* shut down tcl. this also quits the app. */
  SendTCLCommand ( "exit" );

  DebugExitFunction;

  exit( 0 );
}

/* ============================================================ GDF VOLUME */

tkm_tErr LoadGDFHeader ( char* isFSGDHeader,
                         int iRegistrationType,
                         char* isRegistrationName ) {
  MRI* gdInfo = NULL;
  char sTclCommand[tkm_knTclCmdLen];
  char* psTclResult = NULL;
  char sTclResult[tkm_knTclCmdLen];
  int cRead = 0;

  /* Init the GDF display with tcl commands. Save the GDF ID so that
     we can use it to send commands. */
  /*
    set ID [FsgdfPlot_Read $ifnGDF]
    if { $ID < 0 } { return }
    FsgdfPlot_ShowWindow $gGDFID
  */

  sprintf( sTclCommand, "FsgdfPlot_Read %s", isFSGDHeader );
  psTclResult = SendTCLCommand( sTclCommand );
  if ( NULL == psTclResult ) {
    printf( "tkmedit: ERROR: Couldn't read %s\n", isFSGDHeader );
    return tkm_tErr_CouldntLoadGDF;
  }

  strncpy( sTclResult, psTclResult, sizeof(sTclResult) );
  cRead = sscanf( sTclResult, "%d", &gGDFID );
  if ( 1 != cRead ) {
    printf( "tkmedit: ERROR: Couldn't read GDFID from \"%s\"\n", sTclResult );
    return tkm_tErr_CouldntLoadGDF;
  }

  if ( gGDFID < 0 ) {
    printf( "tkmedit: ERROR: FsgdfPlot_Read returned ID %d\n", gGDFID );
    gGDFID = -1;
    return tkm_tErr_CouldntLoadGDF;
  }

  sprintf( sTclCommand, "FsgdfPlot_ShowWindow %d", gGDFID );
  psTclResult = SendTCLCommand( sTclCommand );
  if ( NULL == psTclResult ) {
    return tkm_tErr_CouldntLoadGDF;
  }

  /* now we need to generate a registration transform so get our
     MRIIdx coords into the GDF space. First get the MRI header from
     the GD data file. */
  gdInfo = gdfReadDataInfo( isFSGDHeader );

  /* Now make the registration transform. We'll use this later. */
  MRImakeVox2VoxReg( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues,
                     gdInfo,
                     iRegistrationType,
                     isRegistrationName,
                     &gMRIIdxToGDFIdxTransform );
  if ( NULL == gMRIIdxToGDFIdxTransform ) {
    printf( "tkmedit: ERROR: Couldn't create GDF transform\n" );
    return tkm_tErr_CouldntLoadGDF;
  }

  MRIfree( &gdInfo );

  /* Enable the menu items. */
  tkm_SendTclCommand( tkm_tTclCommand_ShowGDFOptions, "1" );

  return tkm_tErr_NoErr;
}

tkm_tErr GDFPlotMRIIdx ( xVoxelRef iMRIIdx ) {

  tkm_tErr eResult = tkm_tErr_NoErr;

  if ( gGDFID < 0 )
    return tkm_tErr_NoErr;

  /* Just being a list, add this single point, and end the list. */
  eResult = BeginGDFPointList();
  if ( tkm_tErr_NoErr != eResult )
    return eResult;

  eResult = AddGDFPlotMRIIdx( iMRIIdx );
  if ( tkm_tErr_NoErr != eResult )
    return eResult;

  eResult = EndGDFPointList();
  if ( tkm_tErr_NoErr != eResult )
    return eResult;

  return tkm_tErr_NoErr;
}

tkm_tErr BeginGDFPointList () {

  char sTclCommand[tkm_knTclCmdLen];
  char* psTclResult = NULL;

  if ( gGDFID < 0 )
    return tkm_tErr_NoErr;

  /*
    FsgdfPlot_BeginPointList $gGDFID
  */
  sprintf( sTclCommand, "FsgdfPlot_BeginPointList %d", gGDFID );
  psTclResult = SendTCLCommand( sTclCommand );
  if ( NULL == psTclResult ) {
    return tkm_tErr_ErrorAccessingGDF;
  }

  return tkm_tErr_NoErr;
}

tkm_tErr AddGDFPlotMRIIdx ( xVoxelRef iMRIIdx ) {

  char sTclCommand[tkm_knTclCmdLen];
  char* psTclResult = NULL;
  Trns_tErr eTrns = Trns_tErr_NoErr;
  xVoxel gdfIdx;

  if ( gGDFID < 0 )
    return tkm_tErr_NoErr;

  /* Transform the point to a GDF index. */
  eTrns = Trns_ConvertAtoB( gMRIIdxToGDFIdxTransform, iMRIIdx, &gdfIdx );
  if ( Trns_tErr_NoErr != eTrns ) {
    printf( "tkmedit: ERROR: Couldn't convert MRI idx %d, %d, %d\n",
            xVoxl_ExpandInt( iMRIIdx ) );
    return tkm_tErr_ErrorAccessingGDF;
  }

  /*
    FsgdfPlot_AddPoint $gGDFID $x $y $x
  */
  sprintf( sTclCommand, "FsgdfPlot_AddPoint %d %d %d %d",
           gGDFID, xVoxl_ExpandInt( &gdfIdx ) );
  psTclResult = SendTCLCommand( sTclCommand );
  if ( NULL == psTclResult ) {
    return tkm_tErr_ErrorAccessingGDF;
  }

  return tkm_tErr_NoErr;
}

tkm_tErr EndGDFPointList () {

  char sTclCommand[tkm_knTclCmdLen];
  char* psTclResult = NULL;

  if ( gGDFID < 0 )
    return tkm_tErr_NoErr;

  /*
    FsgdfPlot_EndPointList $gGDFID
  */
  sprintf( sTclCommand, "FsgdfPlot_EndPointList %d", gGDFID );
  psTclResult = SendTCLCommand( sTclCommand );
  if ( NULL == psTclResult ) {
    return tkm_tErr_ErrorAccessingGDF;
  }

  return tkm_tErr_NoErr;
}

void tkm_SelectGDFMRIIdx ( xVoxelRef iMRIIdx ) {

  if ( gGDFID < 0 )
    return;

  GDFPlotMRIIdx( iMRIIdx );
}

#ifndef Windows_NT
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
  int tcl_interactive = 0;

  count = read(fileno(stdin), input, BUFFER_SIZE);
  if (count <= 0) {
    if (!gotPartial) {
      if (tty) {
        Tcl_Eval(interp, "exit");
        exit(1);
      } else     {
        Tk_DeleteFileHandler(0);
      }
      return;
    } else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0) {
    if ((input[count-1] != '\n') && (input[count-1] != ';')) {
      gotPartial = 1;
      goto prompt;
    }
    if (!Tcl_CommandComplete(cmd)) {
      gotPartial = 1;
      goto prompt;
    }
  }
  gotPartial = 0;
  Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
  code = Tcl_RecordAndEval(interp, cmd, TCL_EVAL_GLOBAL);
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  Tcl_DStringFree(&command);
  {
    /*
     * Check whether tcl_interactive has been set
     * to enforce interactive-like behavior
     */
    char *ints  = (char *)Tcl_GetVar(interp,
                                     "tcl_interactive",
                                     TCL_GLOBAL_ONLY);
    if (!ints || sscanf(ints,"%d",&tcl_interactive) != 1)
      tcl_interactive = 0;
  }
  if (*interp->result != 0)
    if ((code != TCL_OK) || (tty) || (tcl_interactive) )
      puts(interp->result);
 prompt:
  if (tty) Prompt(interp, gotPartial);
  fflush(stdout);
  Tcl_ResetResult(interp);
}
#endif // Windows_NT

/*=== from TkMain.c ===================================================*/
static void Prompt(interp, partial)
  Tcl_Interp *interp;
  int partial;
{
  char *promptCmd;
  int code;

  promptCmd = (char*) Tcl_GetVar
    (interp,
     partial ? "tcl_prompt2" : "tcl_prompt1",
     TCL_GLOBAL_ONLY);
  if (promptCmd == NULL) {
  defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  } else {
    code = Tcl_Eval(interp, promptCmd);
    if (code != TCL_OK) {
      Tcl_AddErrorInfo(interp,
                       "\n     (script that generates prompt)");
      fprintf(stderr, "%s\n", interp->result);
      goto defaultPrompt;
    }
  }
}

/* ================================================ Control point utilities */

tkm_tErr InitControlPointList () {

  tkm_tErr   eResult = tkm_tErr_NoErr;
  x3Lst_tErr e3DList = xUndL_tErr_NoErr;

  DebugEnterFunction( ("InitControlPointList") );

  /* Free existing list. */
  if ( NULL != gControlPointList ) {
    x3Lst_Delete( &gControlPointList );
  }

  /* Make the list. Conformed volumes are 256^3 so that's our size. */
  e3DList = x3Lst_New( &gControlPointList, kControlPointSpaceDimension );
  DebugAssertThrow( (x3Lst_tErr_NoErr == e3DList) );

  /* Set the comparator for removing duplicates. */
  x3Lst_SetComparator( gControlPointList, CompareVoxels );

  DebugCatch;
  DebugCatchError( e3DList, x3Lst_tErr_NoErr, x3Lst_GetErrorString );
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

void ReadControlPointFile () {

  tkm_tErr  eResult           = tkm_tErr_NoErr;
  Volm_tErr eVolume           = Volm_tErr_NoErr;
  char      sFileName[tkm_knPathLen] = "";
  int       nNumControlPoints = 0;
  int       bUseRealRAS       = 0;
  MPoint*   pControlPoints    = NULL;
  int       nPoint            = 0;
  xVoxel    ras;
  xVoxel    MRIIdx;

  DebugEnterFunction( ("ProcessControlPointFile ()") );

  /* Only load control points for conformed volumes. */
  if ( !mriConformed( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ) ) {
    DebugGotoCleanup;
  }

  /* Make sure our control point list is ready. */
  DebugNote( ("Checking if control point list is inited") );
  DebugAssertThrowX( (NULL != gControlPointList),
                     eResult, tkm_tErr_ErrorAccessingList );

  /* Make the file name */
  DebugNote( ("Making file name from control.dat") );
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints,
                sFileName, sizeof(sFileName) );

  /* If the file doesn't exist, quietly return. */
  DebugNote( ("Checking if %s exists", sFileName) );
  if ( !FileExists( sFileName ) ) {
    DebugQuietThrow();
  }

  /* Read the file. */
  DebugNote( ("Calling MRIreadControlPoints %s", sFileName) );
  pControlPoints =
    MRIreadControlPoints( sFileName, &nNumControlPoints, &bUseRealRAS );
  DebugAssertThrowX( (NULL != pControlPoints), eResult,
                     tkm_tErr_ErrorAccessingFile );

  /* Parse the points. For each one, transform it from real RAS or
     tkReg RAS to MRI index, and add the control point. */
  for ( nPoint = 0; nPoint < nNumControlPoints; nPoint++ ) {

    /* Transform from ras to voxel */
    xVoxl_SetFloat( &ras,
                    pControlPoints[nPoint].x,
                    pControlPoints[nPoint].y,
                    pControlPoints[nPoint].z );

    if ( bUseRealRAS ) {
      DebugNote( ("Converting control point from RAS to MRI Idx") );
      eVolume =
        Volm_ConvertRASToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                 &ras, &MRIIdx );
    } else {
      DebugNote( ("Converting control point from surface RAS to MRI Idx") );
      eVolume =
        Volm_ConvertSurfaceRASToMRIIdx(gAnatomicalVolume[tkm_tVolumeType_Main],
                                       &ras, &MRIIdx );
    }
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                       eResult, tkm_tErr_ErrorAccessingVolume );

    /* Add it to our control points list */
    DebugNote( ("Adding control point") );
    AddMRIIdxControlPoint( &MRIIdx, FALSE );
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  if ( NULL != pControlPoints ) {
    free( pControlPoints );
  }

  DebugExitFunction;
}

void WriteControlPointFile () {

  tkm_tErr   eResult                  = tkm_tErr_NoErr;
  Volm_tErr  eVolume                  = Volm_tErr_NoErr;
  x3Lst_tErr e3DList                  = x3Lst_tErr_NoErr;
  xList_tErr eList                    = xList_tErr_NoErr;
  char       sFileName[tkm_knPathLen] = "";
  int        nPlaneSize               = 0;
  int        nPlane                   = 0;
  xListRef   list                     = NULL;
  xVoxelRef  MRIIdx                   = NULL;
  xVoxel     ras;
  int        nNumControlPointsInPlane = 0;
  int        nNumControlPoints        = 0;
  int        nPoint                   = 0;
  MPoint*    pControlPoints           = NULL;

  DebugEnterFunction( ("WriteControlPointFile()") );

  /* Only do control points for conformed volumes. */
  if ( !mriConformed( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ) ) {
    tkm_DisplayError( "Can't Write Control Points",
                      "Volume is not conformed",
                      "You are editing a non-conformed volume, and tkmedit "
                      "can't use control points unless the volume is "
                      "conformed." );
    DebugGotoCleanup;
  }

  /* Make sure our control point list is ready. */
  DebugNote( ("Checking if control point list is inited") );
  DebugAssertThrowX( (NULL != gControlPointList),
                     eResult, tkm_tErr_ErrorAccessingList );

  /* Make the file name */
  DebugNote( ("Making file name from control.dat") );
  MakeFileName( "control.dat", tkm_tFileName_ControlPoints,
                sFileName, sizeof(sFileName) );

  /* If the file name we got is ./control.dat, it means we couldn't
     find a home directory for the subject and the file will be saved
     locally. Warn the user of this. */
  if ( 0 == strcmp( sFileName, "./control.dat" )) {
    if ( !gbGuessWarningSent ) {
      gbGuessWarningSent = TRUE;
      tkm_DisplayError( "No Home Directory",
                        "Couldn't guess home directory",
                        "You specified the -f file to load a subject and "
                        "I can't guess the home directory from the "
                        "path. Therefore, I can't automatically save "
                        "the control points file in the proper place. "
                        "It is now in the directory from which you started "
                        "tkmedit. When you're done, you will need to move "
                        "the file control.dat to the appropriate place "
                        "in the subject directory, usually subject/tmp/." );
    }
  }

  /* Get the ctrl pts in the list... */
  x3Lst_GetPlaneSize( gControlPointList, &nPlaneSize );
  for ( nPlane = 0; nPlane < nPlaneSize; nPlane++ ) {

    /* Get the list for this x value. */
    DebugNote( ("Getting items in x plane %d\n", nPlane) );
    e3DList = x3Lst_GetItemsInXPlane( gControlPointList, nPlane, &list );
    DebugAssertThrowX( (e3DList == x3Lst_tErr_NoErr),
                       eResult, tkm_tErr_ErrorAccessingList );

    /* Count the control points. */
    DebugNote( ("Getting the number of points in this list") );
    xList_GetCount( list, &nNumControlPointsInPlane );
    nNumControlPoints += nNumControlPointsInPlane;
  }

  /* Allocate an MPoints array. */
  DebugNote( ("Allocating array of size %d for control points",
              nNumControlPointsInPlane));
  pControlPoints = calloc( sizeof(MPoint), nNumControlPoints );
  DebugAssertThrowX( (NULL != pControlPoints), eResult,
                     tkm_tErr_CouldntAllocate );

  /* Get the ctrl pts in the list... */
  nPoint = 0;
  for ( nPlane = 0; nPlane < nPlaneSize; nPlane++ ) {

    /* Get the list for this x value. */
    DebugNote( ("Getting items in x plane %d\n", nPlane) );
    e3DList = x3Lst_GetItemsInXPlane( gControlPointList, nPlane, &list );
    DebugAssertThrowX( (e3DList == x3Lst_tErr_NoErr),
                       eResult, tkm_tErr_ErrorAccessingList );

    /* Traverse the list */
    eList = xList_ResetPosition( list );
    void* pvoid = (void*) &MRIIdx;
    while ( (eList = xList_NextFromPos( list, (void**)pvoid ))
            != xList_tErr_EndOfList ) {

      if ( MRIIdx ) {

        /* Transform to ras space. */
        if ( gbUseRealRAS ) {
          DebugNote( ("Transforming MRI Idx to RAS") );
          eVolume =
            Volm_ConvertMRIIdxToRAS(gAnatomicalVolume[tkm_tVolumeType_Main],
                                    MRIIdx, &ras );
        } else {
          DebugNote( ("Transforming MRI Idx to surface RAS") );
          eVolume =
            Volm_ConvertMRIIdxToSurfaceRAS
            ( gAnatomicalVolume[tkm_tVolumeType_Main],
              MRIIdx,
              &ras );
        }
        DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                           eResult, tkm_tErr_ErrorAccessingVolume );

        /* Copy it to the array. */
        pControlPoints[nPoint].x = xVoxl_GetFloatX( &ras );
        pControlPoints[nPoint].y = xVoxl_GetFloatY( &ras );
        pControlPoints[nPoint].z = xVoxl_GetFloatZ( &ras );
        nPoint++;
      }
    }

    DebugAssertThrowX( (eList == xList_tErr_EndOfList),
                       eResult, tkm_tErr_ErrorAccessingList );
  }

  /* Write the control points file. */
  DebugNote( ("Calling MRIwriteControlPoints %s", sFileName) );
  MRIwriteControlPoints( pControlPoints, nNumControlPoints, gbUseRealRAS,
                         sFileName );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  if ( NULL != pControlPoints ) {
    free( pControlPoints );
  }

  DebugExitFunction;
}

float FindNearestMRIIdxControlPoint ( xVoxelRef        iMRIIdx,
                                      mri_tOrientation inPlane,
                                      xVoxelRef*       opCtrlPt ) {

  tBoolean      bFound           = FALSE;
  float         fDistance        = 0;
  float         fClosestDistance = 0;
  float         fFoundDistance   = 0;
  xListRef      list             = NULL;
  xVoxelRef     MRIIdx           = NULL;
  xVoxelRef     closestMRIIdx    = NULL;
  x3Lst_tErr    e3DList          = x3Lst_tErr_NoErr;
  xList_tErr    eList            = xList_tErr_NoErr;

  DebugEnterFunction( ("FindNearestMRIIdxControlPoint( iMRIIdx=%p, inPlane=%d,"
                       "opCtrlPt=%p )", iMRIIdx, (int)inPlane, opCtrlPt) );

  DebugNote( ("Checking params") );
  DebugAssertThrow( (NULL != iMRIIdx) );
  DebugAssertThrow( (inPlane >= 0 && inPlane < mri_knNumOrientations) );
  DebugAssertThrow( (NULL != opCtrlPt) );

  /* Only do control points for conformed volumes. */
  if ( !mriConformed( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ) ) {
    DebugGotoCleanup;
  }

  /* Skip if we don't have a control point list. */
  DebugAssertThrow( (NULL != gControlPointList ) );

  /* Get the list to search in */
  switch ( inPlane ) {
  case mri_tOrientation_Coronal:
    DebugNote( ("Getting items in z plane %d", xVoxl_GetZ(iMRIIdx)) );
    e3DList =
      x3Lst_GetItemsInZPlane( gControlPointList, xVoxl_GetZ(iMRIIdx), &list );
    break;
  case mri_tOrientation_Horizontal:
    DebugNote( ("Getting items in y plane %d", xVoxl_GetY(iMRIIdx)) );
    e3DList =
      x3Lst_GetItemsInYPlane( gControlPointList, xVoxl_GetY(iMRIIdx), &list );
    break;
  case mri_tOrientation_Sagittal:
    DebugNote( ("Getting items in x plane %d", xVoxl_GetX(iMRIIdx)) );
    e3DList =
      x3Lst_GetItemsInXPlane( gControlPointList, xVoxl_GetX(iMRIIdx), &list );
    break;
  default:
    DebugThrow();
  }
  DebugAssertThrow( (x3Lst_tErr_NoErr == e3DList) );

  /* If we got a list... */
  if ( NULL != list ) {

    /* Start with a large distance */
    fClosestDistance = pow( kControlPointSpaceDimension, 3 );
    bFound = FALSE;

    /* Traverse the list */
    eList = xList_ResetPosition( list );
    void* pvoid = (void*) &MRIIdx;
    while ( (eList = xList_NextFromPos( list, (void**)pvoid ))
            != xList_tErr_EndOfList ) {

      if ( MRIIdx ) {

        /* Get the distance to the clicked voxel... */
        fDistance = sqrt(((xVoxl_GetX(iMRIIdx) - xVoxl_GetX(MRIIdx)) *
                          (xVoxl_GetX(iMRIIdx) - xVoxl_GetX(MRIIdx))) +
                         ((xVoxl_GetY(iMRIIdx) - xVoxl_GetY(MRIIdx)) *
                          (xVoxl_GetY(iMRIIdx) - xVoxl_GetY(MRIIdx))) +
                         ((xVoxl_GetZ(iMRIIdx) - xVoxl_GetZ(MRIIdx)) *
                          (xVoxl_GetZ(iMRIIdx) - xVoxl_GetZ(MRIIdx))) );

        /* If it's less than our max, mark the distance and copy the vox */
        if ( fDistance < fClosestDistance ) {
          fClosestDistance = fDistance;
          closestMRIIdx = MRIIdx;
          bFound = TRUE;
        }
      }
    }

    DebugNote( ("Make sure we're at the end of the list") );
    DebugAssertThrow( (xList_tErr_EndOfList == eList) );

    /* If we found a voxel */
    if ( bFound ) {

      /* Return it. */
      *opCtrlPt = closestMRIIdx;
      fFoundDistance = fClosestDistance;
    }
  }

  DebugCatch;
  DebugCatchError( e3DList, x3Lst_tErr_NoErr, x3Lst_GetErrorString );
  DebugCatchError( eList, xList_tErr_NoErr, xList_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return fFoundDistance;
}

void AddMRIIdxControlPoint ( xVoxelRef iMRIIdx,
                             tBoolean  ibWriteToFile ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  MRIIdx  = NULL;

  DebugEnterFunction( ("AddMRIIdxControlPoint( iMRIIdx=%p,ibWriteToFile=%d )",
                       iMRIIdx, ibWriteToFile) );

  DebugNote( ("Checking params") );
  DebugAssertThrow( (NULL != iMRIIdx) );

  /* Only do control points for conformed volumes. */
  if ( !mriConformed( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ) ) {
    DebugGotoCleanup;
  }

  /* Skip if we don't have a control point list. */
  DebugAssertThrow( (NULL != gControlPointList ) );

  /* Allocate a copy of the voxel */
  DebugNote( ("Allocating a voxel") );
  xVoxl_New( &MRIIdx );
  DebugAssertThrow( (NULL != MRIIdx) );

  xVoxl_Copy( MRIIdx, iMRIIdx );

  /* Add the voxel to the ctrl pt space */
  e3DList = x3Lst_AddItem( gControlPointList, MRIIdx, MRIIdx );
  DebugAssertThrow( (x3Lst_tErr_NoErr == e3DList) );

  /* write it to the control point file. */
  if ( ibWriteToFile )
    WriteControlPointFile();
  else
    gbControlPointsChanged = TRUE ;

//  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeDirty, gbControlPointsChanged ?  "1":"0" );
  DebugCatch;
  DebugCatchError( e3DList, x3Lst_tErr_NoErr, x3Lst_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void DeleteMRIIdxControlPoint ( xVoxelRef iMRIIdx,
                             tBoolean  ibWriteToFile ) {

  x3Lst_tErr e3DList = x3Lst_tErr_NoErr;
  xVoxelRef  MRIIdx  = NULL;

  DebugEnterFunction( ("DeleteMRIIdxControlPoint( iMRIIdx=%p )", iMRIIdx) );

  DebugNote( ("Checking params") );
  DebugAssertThrow( (NULL != iMRIIdx) );

  /* Only do control points for conformed volumes. */
  if ( !mriConformed( gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues ) ) {
    DebugGotoCleanup;
  }

  /* Skip if we don't have a control point list. */
  DebugAssertThrow( (NULL != gControlPointList ) );

  /* Point to the voxel. */
  MRIIdx = iMRIIdx;

  /* Find this voxel in the list and delete it. */
  void* pvoid = (void*) &MRIIdx;
  e3DList = x3Lst_RemoveItem( gControlPointList, MRIIdx, (void**)pvoid );
  DebugAssertThrow( (x3Lst_tErr_NoErr == e3DList) );

  /* delete it */
  xVoxl_Delete( &MRIIdx );

  /* write it to the control point file. */
  if ( ibWriteToFile )
    WriteControlPointFile();
  else
    gbControlPointsChanged = TRUE ;

  DebugCatch;
  DebugCatchError( e3DList, x3Lst_tErr_NoErr, x3Lst_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

// ===========================================================================

// ========================================================= SELECTING REGIONS

tkm_tErr InitSelectionModule () {

  gSelectionVolume           = NULL;
  gSelectionVolumeXDimension = 0;
  gSelectionVolumeYDimension = 0;
  gSelectionVolumeZDimension = 0;
  gSelectionCount            = 0;
  AllocateSelectionVolume();

  return tkm_tErr_NoErr;
}

void DeleteSelectionModule () {

  DebugEnterFunction( ("DeleteSelectionModule()") );

  if ( NULL != gSelectionVolume ) {

    DebugNote( ( "Deleting selection volume.\n" ) );
    free( gSelectionVolume );
  }

  DebugExitFunction;
}

tkm_tErr AllocateSelectionVolume () {

  tkm_tErr   eResult  = tkm_tErr_NoErr;
  Volm_tErr  eVolume  = Volm_tErr_NoErr;
  xList_tErr eList    = xList_tErr_NoErr;

  DebugEnterFunction( ("AllocateSelectionVolume()") );

  DebugAssertQuietThrow( (NULL != gAnatomicalVolume[tkm_tVolumeType_Main]) );

  /* If the volume already exists, delete it. */
  if ( NULL != gSelectionVolume )
    Volm_Delete( &gSelectionVolume );

  /* Clone the anatomical. */
  DebugNote( ("Cloning anatomical") );
  eVolume = Volm_DeepClone( gAnatomicalVolume[tkm_tVolumeType_Main],
                            &gSelectionVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_CouldntAllocate );

  DebugNote( ("Setting selection volume to 0") );
  eVolume = Volm_SetAllValues( gSelectionVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* Initialize the selection list. */
  if ( NULL == gSelectionList ) {
    DebugNote( ("Making selection list") );
    eList = xList_New( &gSelectionList );
    DebugAssertThrowX( (xList_tErr_NoErr == eList),
                       eResult, tkm_tErr_CouldntAllocate );
  }

  /* Set it in the window. */
  DebugNote( ("Setting selection list in window.") );
  MWin_SetSelectionSpace( gMeditWindow, -1, gSelectionVolume );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

void AddVoxelsToSelection ( xVoxelRef iaMRIIdx, int inCount ) {

  int        nVoxel  = 0;
  Volm_tErr  eVolume = Volm_tErr_NoErr;
  float      value   = 0;
  xVoxelRef  newVoxel = NULL;

  if ( NULL == gSelectionVolume ) {
    AllocateSelectionVolume();
  }

  /* For each voxel we got... */
  for ( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {

    /* Set this location in the selection volume to 1 */
    eVolume = Volm_GetValueAtMRIIdx( gSelectionVolume,
                                     &iaMRIIdx[nVoxel], &value );
    if ( fabs( 0.0 - value) < 0.00001) {

      Volm_SetValueAtMRIIdx( gSelectionVolume, &iaMRIIdx[nVoxel], 1.0 );
      /* Add it to the list as well. */
      xVoxl_New( &newVoxel );
      xVoxl_Copy( newVoxel, &iaMRIIdx[nVoxel] );
      xList_InsertItem( gSelectionList, (void*)newVoxel );

      /* Inc our selection count. */
      gSelectionCount++;
    }
  }
}

void RemoveVoxelsFromSelection ( xVoxelRef iaMRIIdx, int inCount ) {

  int        nVoxel     = 0;
  Volm_tErr  eVolume    = Volm_tErr_NoErr;
  float      value      = 0;
  xVoxelRef  delVoxel   = NULL;
  xList_tErr eList      = xList_tErr_NoErr;

  if ( NULL == gSelectionVolume )
    return;
  for ( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {

    /* Set this location in the selection volume to 0 */
    eVolume = Volm_GetValueAtMRIIdx( gSelectionVolume,
                                     &iaMRIIdx[nVoxel], &value );
    if ( fabs(1.0 - value) < 0.00001 ) {
      Volm_SetValueAtMRIIdx( gSelectionVolume, &iaMRIIdx[nVoxel], 0 );

      /* Find and remove it from the list. */
      eList = xList_ResetPosition( gSelectionList );
      void* pvoid = (void*) &delVoxel;
      while ( xList_tErr_NoErr ==
              (eList = xList_NextFromPos( gSelectionList, 
                                          (void**)pvoid ))) {
	  if ( NULL != delVoxel ) {
          if ( xVoxl_IsEqualInt( delVoxel, &iaMRIIdx[nVoxel] ) ) {
            void* pvoidVoxel = (void*) &delVoxel;
            xList_RemoveItem( gSelectionList, (void**)pvoidVoxel );
            xVoxl_Delete( &delVoxel );

            /* Dec our selection count. */
            gSelectionCount--;

            break;
          }
        }
      }

    }
  }
}

void ClearSelection () {

  xVoxelRef  delVoxel = NULL;
  xList_tErr eList    = xList_tErr_NoErr;

  if ( NULL == gSelectionVolume )
    return;

  DebugNote( ("Setting selection volume to 0") );
  Volm_SetAllValues( gSelectionVolume, 0 );

  /* Clear the list */
  DebugNote( ("Clearing selection list") );
  void* pvoid = (void*) &delVoxel;
  while ( xList_tErr_NoErr ==
          (eList = xList_PopItem( gSelectionList, (void**)pvoid ))) {
    if ( NULL != delVoxel ) {
      xVoxl_Delete( &delVoxel );
    }
  }

  /* Zero the selection count */
  gSelectionCount = 0;
}

void SaveSelectionToLabelFile ( char * isFileName ) {

  tkm_tErr      eResult                  = tkm_tErr_NoErr;
  Volm_tErr     eVolume                  = Volm_tErr_NoErr;
  char          sFileName[tkm_knPathLen] = "";
  float         value                    = 0;
  xVoxelRef     MRIIdx                   = NULL;
  xList_tErr    eList                    = xList_tErr_NoErr;
  xVoxel        ras;
  LABEL*        pLabel                   = NULL;
  LABEL_VERTEX* pVertex                  = NULL;
  int           nVoxel                   = 0;
  int           eLabel                   = 0;

  DebugEnterFunction( ("SaveSelectionToLabelFile ( isFileName=%s )",
                       isFileName) );

  DebugAssertThrowX( (NULL != gSelectionVolume),
                     eResult, tkm_tErr_NoErr );

  /* make the file name */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_Label,
                sFileName, sizeof(sFileName) );

  /* allocate a label file with that number of voxels, our subject name,
     and the passed in label file name. */
  DebugNote( ("Allocating label with %d voxels") );
  pLabel = LabelAlloc( gSelectionCount, NULL, sFileName );
  DebugAssertThrowX( (NULL != pLabel), eResult, tkm_tErr_CouldntAllocate );

  /* set the number of points in the label */
  pLabel->n_points = gSelectionCount;

  /* Copy subject name. */
  Volm_CopySubjectName( gAnatomicalVolume[tkm_tVolumeType_Main],
                        pLabel->subject_name, 100 );

  nVoxel = 0;

  /* Go through the selection list. */
  eList = xList_ResetPosition( gSelectionList );
  void* pvoid = (void*) &MRIIdx;
  while ( xList_tErr_NoErr ==
          (eList = xList_NextFromPos( gSelectionList, (void**)pvoid ))) {

    if ( NULL != MRIIdx ) {

      /* convert mri idx to surface ras. note we may use surface ras
         here because it ignores c_ras, which is what label files
         should to surfacebe comptaible with tksurfer.  */
      if ( gbUseRealRAS ) {
        eVolume =
          Volm_ConvertMRIIdxToRAS( gSelectionVolume, MRIIdx, &ras );
      } else {
        eVolume =
          Volm_ConvertMRIIdxToSurfaceRAS( gSelectionVolume, MRIIdx, &ras );
      }
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

      /* set the vno to -1, which is significant somewhere outside
         the realm of tkmedit. set stat value to the mri value
         and deleted to not */
      Volm_GetValueAtMRIIdx_( gAnatomicalVolume[tkm_tVolumeType_Main],
                              MRIIdx, &value );
      pVertex->vno = -1;
      pVertex->stat = value;
      pVertex->deleted = FALSE;

      /* inc our global count. */
      nVoxel++;
    }
  }

  /* write the file */
  DebugNote( ("Writing label file to %s", sFileName) );
  eLabel = LabelWrite( pLabel, sFileName );
  DebugAssertThrowX( (NO_ERROR == eLabel),
                     eResult, tkm_tErr_CouldntWriteFile );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  if ( NULL != pLabel ) {
    DebugNote( ("Deleting label") );
    LabelFree( &pLabel );
  }

  DebugExitFunction;
}

tkm_tErr LoadSelectionFromLabelFile ( char* isFileName ) {

  tkm_tErr      eResult         = tkm_tErr_NoErr;
  Volm_tErr     eVolume         = Volm_tErr_NoErr;
  char          sFileName[tkm_knPathLen]   = "";
  LABEL*        pLabel         = NULL;
  LABEL_VERTEX* pVertex         = NULL;
  int           nNumVoxels       = 0;
  int           nVoxel         = 0;
  xVoxel        ras;
  xVoxel        MRIIdx;
  char          sError[tkm_knErrStringLen] = "";

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
  for ( nVoxel = 0; nVoxel < nNumVoxels; nVoxel++ ) {

    /* get the vertex. */
    pVertex = &(pLabel->lv[nVoxel]);

    /* only process verticies that arn't deleted. */
    if ( !(pVertex->deleted) ) {

      /* transform from ras to voxel. note we may use surface ras
         here, see note in SaveSelectionToLabelFile. */
      xVoxl_SetFloat( &ras, pVertex->x, pVertex->y, pVertex->z );
      if ( gbUseRealRAS ) {
        eVolume =
          Volm_ConvertRASToMRIIdx( gSelectionVolume, &ras, &MRIIdx );
      } else {
        eVolume =
          Volm_ConvertSurfaceRASToMRIIdx( gSelectionVolume, &ras, &MRIIdx );
      }
      DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                         eResult, tkm_tErr_ErrorAccessingVolume );

      /* add to selection */
      AddVoxelsToSelection( &MRIIdx, 1 );
    }
  }

  /* set the window's selection again to force a redraw. */
  MWin_SetSelectionSpace( gMeditWindow, -1, gSelectionVolume );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  xUtil_snprintf( sError, sizeof(sError), "Loading label %s", isFileName );
  tkm_DisplayError( sError,
                    tkm_GetErrorString(eResult),
                    "Tkmedit couldn't read the label you "
                    "specified. This could be because the format "
                    "wasn't valid or the file wasn't found." );
  EndDebugCatch;

  if ( NULL != pLabel ) {
    DebugNote( ("Deleting label") );
    LabelFree( &pLabel );
  }

  DebugExitFunction;

  return eResult;
}

void GraphSelectedRegion () {

  FunV_tErr  eFunctional = FunV_tErr_NoError;
  int        nDimensionX = 0;
  int        nDimensionY = 0;
  int        nDimensionZ = 0;
  float      value       = 0;
  xVoxel     idx;
  xVoxel     MRIIdx;

  DebugEnterFunction( ("GraphSelectedRegion()") );

  DebugAssertThrow( (NULL != gSelectionVolume ) );

  // clear the functional display list.
  eFunctional = FunV_BeginSelectionRange( gFunctionalVolume );
  if ( FunV_tErr_NoError != eFunctional )
    goto error;

  /* Look for selected voxels. */
  Volm_GetDimensions( gSelectionVolume, &nDimensionX,
                      &nDimensionY, &nDimensionZ );
  xVoxl_Set( &idx, 0, 0, 0 );
  while ( xVoxl_IncrementUntilLimits( &idx, nDimensionX-1,
                                      nDimensionY-1, nDimensionZ-1 )) {

    Volm_GetValueAtIdxUnsafe( gSelectionVolume, &idx, &value );
    if ( 0 != value ) {

      /* add it to the functional display list. */
      Volm_ConvertIdxToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                               &idx, &MRIIdx );
      eFunctional =
        FunV_AddMRIIdxToSelectionRange( gFunctionalVolume, &MRIIdx );
      if ( FunV_tErr_NoError != eFunctional )
        goto error;
    }
  }

  /* finish the list */
  eFunctional = FunV_EndSelectionRange( gFunctionalVolume );
  if ( FunV_tErr_NoError != eFunctional )
    goto error;

  DebugCatch;

  if ( eFunctional != FunV_tErr_NoError ) {
    DebugPrint( ( "FunV error %d in GraphSelectedRegion.\n", eFunctional ) );
    OutputPrint "Error graphing selection.\n" EndOutputPrint;
  }

  EndDebugCatch;

  DebugExitFunction;
}

void SelectVoxelsByFuncValue ( FunV_tFindStatsComp iCompare ) {

  xVoxel cursor;

  DebugEnterFunction( ("SelectVoxelsByFuncValue") );

  if ( NULL == gMeditWindow ||
       NULL == gFunctionalVolume )
    return;

  /* get cursor */
  DebugNote( ("Getting cursor") );
  MWin_GetCursorInMRIIdx( gMeditWindow, &cursor );

  /* select voxels */
  DebugNote( ("Selecting voxels around %d,%d,%d", xVoxl_ExpandInt( &cursor )));
  FunV_FloodSelect( gFunctionalVolume, &cursor, tkm_tVolumeType_Main,
                    9999, iCompare );

  DebugExitFunction;
}

tkm_tErr FloodSelect ( xVoxelRef         iSeedAnaIdx,
                       tBoolean          ib3D,
                       tkm_tVolumeTarget iSrc,
                       float             iFuzzy,
                       float             iDistance,
                       tBoolean          ibSelect ) {

  tkm_tErr                    eResult      = tkm_tErr_NoErr;
  Volm_tFloodParams           params;
  tkm_tFloodSelectCallbackData callbackData;
  mriVolumeRef                sourceVolume = NULL;
  Volm_tErr                   eVolume      = Volm_tErr_NoErr;
  char                        sTclArguments[tkm_knTclCmdLen] = "";
  xVoxel                      mriIdx;

  DebugEnterFunction( ("FloodSelect( iSeedAnaIdx=%p, ib3D=%d, iSrc=%d, "
                       "iFuzzy=%f, iDistance=%f", iSeedAnaIdx, ib3D,
                       iSrc, iFuzzy, iDistance) );

  DebugAssertThrowX( (NULL != iSeedAnaIdx),
                     eResult, tkm_tErr_InvalidParameter );

  /* We need an MRI index. */
  Volm_ConvertIdxToMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                           iSeedAnaIdx, &mriIdx );

  xVoxl_Copy( &params.mSourceIdx, &mriIdx );
  params.mfFuzziness             = iFuzzy;
  params.mComparatorType         = Volm_tValueComparator_EQ;
  params.mComparatorFunc         = NULL;
  params.mfMaxDistance           = iDistance;
  params.mb3D                    = ib3D;
  MWin_GetOrientation ( gMeditWindow, &params.mOrientation );

  /* Set the callback function data. Tell it to use the callback data
     we just initialized. */
  params.mpFunction     = FloodSelectCallback;
  params.mpFunctionData = (void*)&callbackData;

  callbackData.mbSelect = ibSelect;
  callbackData.mnCount  = 0;

  /* See what source volume to use. For everything other than the main
     anatomical volume, check to see if it exists first. For the
     segmentation targets, set the fuzzinees to 0, since it doesn't
     really apply to label values. */
  switch ( iSrc ) {
  case tkm_tVolumeTarget_MainAna:
    sourceVolume = gAnatomicalVolume[tkm_tVolumeType_Main];
    break;
  case tkm_tVolumeTarget_AuxAna:
    if ( NULL == gAnatomicalVolume[tkm_tVolumeType_Aux] ) {
      strcpy( sTclArguments,
              "\"Cannot use aux volume as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gAnatomicalVolume[tkm_tVolumeType_Aux];
    break;
  case tkm_tVolumeTarget_MainSeg:
    if ( NULL == gSegmentationVolume[tkm_tSegType_Main] ) {
      strcpy( sTclArguments,
              "\"Cannot use segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gSegmentationVolume[tkm_tSegType_Main];
    params.mfFuzziness = 0;
    break;
  case tkm_tVolumeTarget_AuxSeg:
    if ( NULL == gSegmentationVolume[tkm_tSegType_Aux] ) {
      strcpy
        ( sTclArguments,
          "\"Cannot use aux segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gSegmentationVolume[tkm_tSegType_Aux];
    params.mfFuzziness = 0;
    break;
  default:
    DebugThrowX( eResult, tkm_tErr_InvalidParameter );
    break;
  }

  /* Now get the source volume value. */
  Volm_GetValueAtMRIIdx( sourceVolume, &mriIdx, &params.mfSourceValue );

  /* Start listening for a cancel. */
  xUtil_StartListeningForUserCancel();

  /* Do it! */
  eVolume = Volm_Flood( sourceVolume, &params );

  /* If we selected more than 1000 voxels, we printed a message and
     started printing update dots. Now close off the message. */
  if ( callbackData.mnCount > 1000 ) {
    printf( "done. %d voxels selected. \n", callbackData.mnCount );
  }

  /* Stop listening for the cancel. */
  xUtil_StopListeningForUserCancel();

  UpdateAndRedraw();

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tVisitCommand FloodSelectCallback ( xVoxelRef iMRIIdx,
                                         float     iValue,
                                         void*     iData ) {

  tkm_tFloodSelectCallbackData* callbackData;

  callbackData = (tkm_tFloodSelectCallbackData*)iData;

  /* Incremenet our count. If it's over 1000, print a message saying
     the user can cancel and start printing update dots. */
  callbackData->mnCount++;
  if ( callbackData->mnCount == 1000 ) {
    printf( "Selecting (press ctrl-c to cancel) " );
  }
  if ( callbackData->mnCount > 1000 &&
       callbackData->mnCount % 100 == 0 ) {
    printf( "." );
    fflush( stdout );
  }

  /* Check the user cancel. If they canceled, stop. */
  if ( xUtil_DidUserCancel() ) {
    return Volm_tVisitComm_Stop;
  }

  /* Select or deselect this voxel. */
  if ( callbackData->mbSelect ) {
    AddVoxelsToSelection( iMRIIdx, 1 );
  } else {
    RemoveVoxelsFromSelection( iMRIIdx, 1 );
  }

  return Volm_tVisitComm_Continue;
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
  tkm_SendTclCommand
    ( tkm_tTclCommand_UpdateSubjectDirectory, gsSubjectHomeDir );

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
  if ( NULL == isInput )
    return;

  /* look at the first character */
  switch ( isInput[0] ) {

    /* tilde, attach the subject home dir. */
  case '~':

    /* look at the second char. if it's a / or null use this subject. */
    if ( isInput[1] == '/' ||
         isInput[1] == '\0' ) {

      xUtil_snprintf( sFileName, sizeof(sFileName),
                      "%s%s", gsSubjectHomeDir, &(isInput[1]) );

      /* if it's not, treat the part after it as the subject name */
    } else {
      xUtil_snprintf( sFileName, sizeof(sFileName),
                      "%s/../%s", gsSubjectHomeDir, &(isInput[1]) );
    }

    break;

    /* dot. could be ./ or ../ . if the former, replace the dot with
       the cwd, if the latter, prepend the cwd and a slash. */
  case '.':

    if ( isInput[1] == '.' ) {
      /* ../ -> cwd/filename */
      xUtil_snprintf( sFileName, sizeof(sFileName),
                      "%s/%s", gsUserHomeDir, isInput );
    } else {
      /* ./ -> cwdfilename[1] */
      xUtil_snprintf( sFileName, sizeof(sFileName),
                      "%s%s", gsUserHomeDir, &(isInput[1]) );
    }
    break;

    /* slash, it's a full path. */
  case '/':
    xUtil_strncpy( sFileName, isInput, sizeof(sFileName) );
    break;

    /* else, prepend subject home, then sub dir, then rest of file name. */
  default:
    if ( gEnableFileNameGuessing ) {
      if ( iType != tkm_tFileName_PWD ) {
        xUtil_snprintf( sFileName, sizeof(sFileName),
                        "%s/%s/%s", gsSubjectHomeDir,
                        ksaFileNameSubDirs[iType], isInput );
      } else {
        xUtil_snprintf(sFileName,sizeof(sFileName), "./%s", isInput);
      }
    } else {
      xUtil_snprintf(sFileName,sizeof(sFileName), "./%s", isInput);
      /* xUtil_strncpy( sFileName, isInput, sizeof(sFileName) ); */
    }
    break;
  }

  xUtil_strncpy( osCompleteFileName, sFileName, inDestSize );

  DebugExitFunction;
}

// ===========================================================================


/* ===================================================== General utilities */

void SetTclInterp ( Tcl_Interp * inInterp ) {

  gTclInterp = inInterp;
}

Tcl_Interp * GetTclInterp () {

  return gTclInterp;
}

char* SendTCLCommand ( char * inCommand ) {

  int rTcl;
  Tcl_Interp * theInterp;
  char* sTclResult = NULL;
  char sTclCmd[tkm_knTclCmdLen];
  xGArr_tErr eList = xGArr_tErr_NoErr;

  // get the interp and send the command.
  theInterp = GetTclInterp ();

  if ( NULL != theInterp ) {

    rTcl = Tcl_Eval( theInterp, inCommand );
    sTclResult = (char *)Tcl_GetStringResult( theInterp );

    // print any error msgs.
    if ( TCL_OK != rTcl ) {
      DebugPrint( ( "Cmd: %s\n", inCommand ) );
      DebugPrint( ( "\tCommand did not return OK.\n" ) );
      DebugPrint( ( "\tResult: %s\n", theInterp->result ) );
    }

  } else {

    if ( NULL == gCachedTclCommands ) {
      eList = xGArr_New( &gCachedTclCommands, sizeof( sTclCmd ), 20 );
      if ( xGArr_tErr_NoErr != eList ) {
        DebugPrint( ( "Couldn't allocate list: couldn't cache %s\n",
                      inCommand ) );
        return NULL;
      }
    }

    /* cache the command so we can send it later */
    strcpy( sTclCmd, inCommand );
    eList = xGArr_Add( gCachedTclCommands, sTclCmd );
    if ( xGArr_tErr_NoErr != eList ) {
      DebugPrint( ( "Couldn't cache %s\n", inCommand ) );
    }
  }

  return sTclResult;
}

void SendCachedTclCommands () {

  char sTclCmd[tkm_knTclCmdLen];
  xGArr_tErr eList = xGArr_tErr_NoErr;

  if ( NULL == gCachedTclCommands )
    return;

  /* send all our cached commands */
  eList = xGArr_ResetIterator( gCachedTclCommands );
  if ( xGArr_tErr_NoErr != eList )
    goto error;
  while ( (eList = xGArr_NextItem( gCachedTclCommands, (void*)&sTclCmd ))
          == xGArr_tErr_NoErr ) {

    SendTCLCommand( sTclCmd );
  }

  goto cleanup;
 error:

  /* print error message */
  if ( xGArr_tErr_NoErr != eList ) {
    DebugPrint( ("Error %d in SendTCLCommand: %s\n",
                 eList, xGArr_GetErrorString(eList) ) );
  }
 cleanup:

  xGArr_Delete( &gCachedTclCommands );
}


// ============================================================= VOLUME ACCESS


tkm_tErr LoadVolume ( tkm_tVolumeType iType,
                      char*           isName,
                      tBoolean        ibConform ) {

  tkm_tErr             eResult                        = tkm_tErr_NoErr;
  Volm_tErr            eVolume                        = Volm_tErr_NoErr;
  MWin_tErr            eWindow                        = MWin_tErr_NoErr;
  char                 sPath[tkm_knPathLen]           = "";
  char*                pEnd                           = NULL;
  char                 sError[tkm_knErrStringLen]     = "";
  mriVolumeRef         newVolume                      = NULL;
  char                 sTclArguments[tkm_knTclCmdLen] = "";
  Volm_tSampleType     sampleType           = Volm_tSampleType_Nearest;
  Volm_tResampleMethod resampleMethod       = Volm_tResampleMethod_RAS;
  float                fBrightness          = Volm_kfDefaultBrightness;
  float                fContrast            = Volm_kfDefaultContrast;
  int                  mainDimensions[3]    = {0, 0, 0};
  int                  auxDimensions[3]     = {0, 0, 0};

  DebugEnterFunction( ("LoadVolume( iType=%d, isName=%s, ibConform = %d )",
                       (int)iType, isName, ibConform) );

  DebugAssertThrowX( (NULL != isName), eResult, tkm_tErr_InvalidParameter );

  /* make a filename */
  DebugNote( ("Making filename from %s", isName) );
  MakeFileName( isName, tkm_tFileName_Volume, sPath, sizeof(sPath) );

  /* if there is a /COR- at the end, take it off. */
  DebugNote( ("Looking for COR- at end of file name") );
  pEnd = strstr( sPath, "/COR-" );
  if ( NULL != pEnd )
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

  /* Conform if desired. */
  if ( ibConform ) {
    Volm_Conform( newVolume );
  }

  /* If this is the aux volume, make sure it is the same size as the
     main volume. */
  if ( tkm_tVolumeType_Aux == iType ) {
    Volm_GetDimensions( gAnatomicalVolume[tkm_tVolumeType_Main],
                        &mainDimensions[0], &mainDimensions[1],
                        &mainDimensions[2] );
    Volm_GetDimensions( newVolume,
                        &auxDimensions[0], &auxDimensions[1],
                        &auxDimensions[2] );
    DebugAssertThrowX( (mainDimensions[0] == auxDimensions[0] &&
                        mainDimensions[1] == auxDimensions[1] &&
                        mainDimensions[2] == auxDimensions[2]),
                       eResult, tkm_tErr_MainAuxVolumesDifferentSize );
  }

  if ( gbScaleUpVolume ) {
    Volm_SetMinVoxelSizeToOne( newVolume );
    Volm_GetMRIIdxToAnaIdxTransform( newVolume,
                                     &gMRIIdxToAnaIdxTransform );
  }

  /* if the volume exists, get the brightness and contrast to restore
     later, and delete the volume */
  if ( NULL != gAnatomicalVolume[iType] ) {
    Volm_GetBrightnessAndContrast( gAnatomicalVolume[iType],
                                   &fBrightness, &fContrast );
    UnloadDisplayTransform( iType );
    Volm_Delete( &gAnatomicalVolume[iType] );
  }

  /* save the new volume */
  gAnatomicalVolume[iType] = newVolume;

  /* show the tal coords and hide the ras coords */
  if (NULL != gAnatomicalVolume[iType]->mpMriValues->linear_transform) {
    DebugNote( ("Showing Tal coords") );
    tkm_SendTclCommand( tkm_tTclCommand_ShowTalCoords, "1" );
    DebugNote( ("Hiding RAS coords") );
    tkm_SendTclCommand( tkm_tTclCommand_ShowRASCoords, "0" );
  } else {
    DebugNote( ("Hiding Tal coords") );
    tkm_SendTclCommand( tkm_tTclCommand_ShowTalCoords, "0" );
    DebugNote( ("Showing RAS coords") );
    tkm_SendTclCommand( tkm_tTclCommand_ShowRASCoords, "1" );
  }

  /* AnaIdx coordinate dimensions. */
  gnAnatomicalDimensionX = 256;
  gnAnatomicalDimensionY = 256;
  gnAnatomicalDimensionZ = 256 ;

  /* Set the default color scale. Get the value min and max from the
     volume and use that. Set brightness and contrast, which will be
     the default values, or, if we had a volume loaded before, the
     values from that volume. */
  Volm_GetValueMinMax( gAnatomicalVolume[iType],
                       &gfaAnaColorMin[iType], &gfaAnaColorMax[iType] );
  sprintf( sTclArguments, "%d %.2f %.2f", (int)iType,
           gfaAnaColorMin[iType], gfaAnaColorMax[iType] );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValueMinMax, sTclArguments );


  SetVolumeColorMinMax( iType, gfaAnaColorMin[iType], gfaAnaColorMax[iType] );
  SetVolumeBrightnessAndContrast( iType, fBrightness, fContrast );

  /* volume is clean */
  gbAnatomicalVolumeDirty[iType] = FALSE;

  /* If this is the main volume, we need to initialize new control
     point and selection stuff, because they are based on the
     dimensions of the main volume. Also set them in the window. */
  if ( tkm_tVolumeType_Main == iType ) {
    DebugNote( ("Initializing control point list") );
    eResult = InitControlPointList();
    DebugAssertThrow( (eResult == tkm_tErr_NoErr) );

    /* Reread the control points. */
    ReadControlPointFile();

    DebugNote( ("Setting control points space in window") );
    eWindow =
      MWin_SetControlPointsSpace( gMeditWindow, -1, gControlPointList );
    DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );

    DebugNote( ("Allocating selection volume") );
    eResult = AllocateSelectionVolume();
    DebugAssertThrow( (eResult == tkm_tErr_NoErr) );

    DebugNote( ("Setting selection list in window.") );
    eWindow = MWin_SetSelectionSpace( gMeditWindow, -1, gSelectionVolume );
    DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );
  }

  /* set data in window */
  if ( NULL != gMeditWindow ) {
    DebugNote( ("Setting volume in main window") );
    switch ( iType ) {
    case tkm_tVolumeType_Main:
      eWindow = MWin_SetVolume( gMeditWindow, -1,
                                gAnatomicalVolume[iType],
                                gnAnatomicalDimensionX,
                                gnAnatomicalDimensionY,
                                gnAnatomicalDimensionZ  );
      DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );

      /* get a ptr to the idx to ras transform */
      DebugNote( ("Getting a pointer to the MRI idx to ana idx transform") );
      Volm_GetMRIIdxToAnaIdxTransform( gAnatomicalVolume[iType],
                                       &gMRIIdxToAnaIdxTransform );
      break;
    case tkm_tVolumeType_Aux:
      eWindow =  MWin_SetAuxVolume( gMeditWindow, -1,
                                    gAnatomicalVolume[iType],
                                    gnAnatomicalDimensionX,
                                    gnAnatomicalDimensionY,
                                    gnAnatomicalDimensionZ );
      DebugAssertThrow( (MWin_tErr_NoErr == eWindow ) );

      break;
    default:
      break;
    }
  }

  gm_screen2ras = extract_i_to_r( gAnatomicalVolume[iType]->mpMriValues );
  gm_ras2screen = extract_r_to_i( gAnatomicalVolume[iType]->mpMriValues );

  /* If this is the main volume and it's a COR, turn on the
     axes. Otherwise turn them off. */
  if ( iType == tkm_tVolumeType_Main &&
       gAnatomicalVolume[iType]->mpMriValues->x_r == -1.0 &&
       gAnatomicalVolume[iType]->mpMriValues->x_a ==  0.0 &&
       gAnatomicalVolume[iType]->mpMriValues->x_s ==  0.0 &&
       gAnatomicalVolume[iType]->mpMriValues->y_r ==  0.0 &&
       gAnatomicalVolume[iType]->mpMriValues->y_a ==  0.0 &&
       gAnatomicalVolume[iType]->mpMriValues->y_s == -1.0 &&
       gAnatomicalVolume[iType]->mpMriValues->z_r ==  0.0 &&
       gAnatomicalVolume[iType]->mpMriValues->z_a ==  1.0 &&
       gAnatomicalVolume[iType]->mpMriValues->z_s ==  0.0 ) {
    MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_Axes, TRUE );
  } else {
    MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_Axes, FALSE );
  }

  /* Send info to the tcl window. */
  Volm_GetResampleMethod( gAnatomicalVolume[iType], &resampleMethod );
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %d",
                  (int)iType, (int)resampleMethod );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeResampleMethod,
                      sTclArguments );

  Volm_GetSampleType( gAnatomicalVolume[iType], &sampleType );
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %d",
                  (int)iType, (int)sampleType );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeSampleType, sTclArguments );

  /* Enable control points only if volume is conformed. If the user
     overrided this, enable them. */
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d",
                  mriConformed( gAnatomicalVolume[iType]->mpMriValues ) );
  if( gbForceEnableControlPoints ) {
    xUtil_strncpy( sTclArguments, "1", sizeof(sTclArguments) );
  }
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeIsConformed,
                      sTclArguments );
  tkm_SendTclCommand( tkm_tTclCommand_ShowControlPointsOptions,
                      sTclArguments );


  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  xUtil_snprintf( sError, sizeof(sError), "Loading volume %s", isName );
  tkm_DisplayError( sError, tkm_GetErrorString(eResult),
                    "Tkmedit couldn't read the volume you specified.\n"
                    "  This could be because the image format wasn't "
                    "recognized,\n  or it couldn't find the proper header,\n"
                    "  or the file(s) were unreadable,\n  or it was the wrong "
                    "size." );

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


  switch ( iType ) {
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
  if ( NULL != gAnatomicalVolume[iType] ) {
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
  if ( NULL != gAnatomicalVolume[iVolume] ) {
    eVolume = Volm_LoadDisplayTransform( gAnatomicalVolume[iVolume],
                                         sFileName );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume ),
                       eResult, tkm_tErr_ErrorAccessingVolume );
  }

  /* enable the menu options */
  switch ( iVolume ) {
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
  if ( NULL != gMeditWindow )
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
  if ( NULL != gAnatomicalVolume[iVolume] ) {
    eVolume = Volm_UnloadDisplayTransform( gAnatomicalVolume[iVolume] );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume ),
                       eResult, tkm_tErr_ErrorAccessingVolume );
  }

  /* disable the menu options */
  switch ( iVolume ) {
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
  if ( NULL != gMeditWindow )
    MWin_RedrawAll( gMeditWindow );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

tkm_tErr SaveVolume ( tkm_tVolumeType iVolume,
                      char*           isPath,
                      tBoolean        ibSetFileName) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  tBoolean  bDirty         = FALSE;
  char      sFileName[tkm_knPathLen] = "";
  char*     psFileName         = NULL;
  char      sBaseFileName[tkm_knPathLen] = "";
  char      sTouchFileName[tkm_knPathLen] = "";
  FILE*     fTouch = NULL;

  DebugEnterFunction( ("SaveVolume( iVolume=%d, isPath=%s )",
                       iVolume, isPath) );

  if ( NULL == gAnatomicalVolume[iVolume] )
    DebugGotoCleanup;

  /* if we have editing disabled, return */
  if ( !editflag ) {
    tkm_DisplayError( "Saving Volume",
                      "Couldn't save volume",
                      "This session was started in read-only mode. You cannot "
                      "save changes." );
    DebugGotoCleanup;
  }

  /* if we haven't been edited, return */
  eResult = IsVolumeDirty( iVolume, &bDirty );
  if ( !bDirty ) {
    DebugPrint( ("SaveVolume called when volume is clean!\n") );
    DebugGotoCleanup;
  }

  /* make a file name */
  if ( NULL != isPath ) {
    DebugNote( ("Making file name from %s", isPath) );
    MakeFileName( isPath, tkm_tFileName_Volume, sFileName, sizeof(sFileName) );
    psFileName = sFileName;
  }

  /* write the main anatomical volume */
  eVolume = Volm_Save( gAnatomicalVolume[iVolume], psFileName, ibSetFileName );
  if ( Volm_tErr_NoErr != eVolume ) {
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

  /* Try to create the touch file. If ../touch doesn't exist, this
     will fail, but that's okay because if ../touch doesn't exist, we
     don't need to create the file. */
  Volm_CopyVolumeName( gAnatomicalVolume[iVolume],
                       sBaseFileName, sizeof(sBaseFileName) );
  MakeFileName( sBaseFileName, tkm_tFileName_Touch,
                sTouchFileName, sizeof(sTouchFileName) );
  sprintf( sTouchFileName, "%s.tkmedit.touch", sTouchFileName );

  fTouch = fopen( sTouchFileName, "w" );
  if ( fTouch ) {
    fclose( fTouch );
  }

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
  switch ( iVolume ) {
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
    break;
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
                                      float           ifBrightness,
                                      float           ifContrast ) {

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
  SendVolumeColorScaleUpdate( iVolume );

  /* big redraw */
  MWin_RedrawAll( gMeditWindow );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void SetVolumeColorMinMax ( tkm_tVolumeType iVolume,
                            float           ifMin,
                            float           ifMax ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;

  DebugEnterFunction( ("SetVolumeColorMinMax ( iVolume=%d, "
                       "ifMin=%.2f, ifMax=%.2f )",
                       (int)iVolume, ifMin, ifMax) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* save these values. */
  gfaAnaColorMin[iVolume] = ifMin;
  gfaAnaColorMax[iVolume] = ifMax;

  /* set the values in the volume. */
  DebugNote( ("Setting min and max") );
  eVolume = Volm_SetColorMinMax( gAnatomicalVolume[iVolume],
                                 gfaAnaColorMin[iVolume],
                                 gfaAnaColorMax[iVolume] );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* update the tcl window */
  SendVolumeColorScaleUpdate( iVolume );

  /* big redraw */
  MWin_RedrawAll( gMeditWindow );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void SetVolumeSampleType  ( tkm_tVolumeType  iVolume,
                            Volm_tSampleType iType ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("SetVolumeSampleType ( iVolume=%d, iTyoe=%d )",
                       (int)iVolume, (int)iType) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* Set the type in the volume. */
  DebugNote( ("Setting sample type") );
  eVolume = Volm_SetSampleType( gAnatomicalVolume[iVolume], iType );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* update the tcl window */
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %d",
                  (int)iVolume, (int)iType );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeSampleType, sTclArguments );

  /* big redraw */
  MWin_RedrawAll( gMeditWindow );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void SetVolumeResampleMethod  ( tkm_tVolumeType      iVolume,
                                Volm_tResampleMethod iMethod ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("SetVolumeResampleMethod ( iVolume=%d, iMethod=%d )",
                       (int)iVolume, (int)iMethod) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* Set the type in the volume. */
  DebugNote( ("Setting resample method") );
  eVolume = Volm_SetResampleMethod( gAnatomicalVolume[iVolume], iMethod );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* If this is the main volume, also set the selection volume. */
  if ( tkm_tVolumeType_Main == iVolume ) {
    DebugNote( ("Setting resample method in selection volume") );
    eVolume = Volm_SetResampleMethod( gSelectionVolume, iMethod );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                       eResult, tkm_tErr_ErrorAccessingVolume );
  }

  /* update the tcl window */
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %d",
                  (int)iVolume, (int)iMethod );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeResampleMethod,
                      sTclArguments );

  /* big redraw */
  MWin_RedrawAll( gMeditWindow );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void SendVolumeColorScaleUpdate ( tkm_tVolumeType iVolume ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  char      sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("SendVolumeColorScaleUpdate ( iVolume=%d )",
                       (int)iVolume) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* update the tcl window */
  DebugNote( ("Sending color scale update to tcl window") );
  xUtil_snprintf( sTclArguments, sizeof(sTclArguments), "%d %f %f %f %f",
                  (int)iVolume, gfaBrightness[iVolume], gfaContrast[iVolume],
                  gfaAnaColorMin[iVolume], gfaAnaColorMax[iVolume]);
  tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeColorScale, sTclArguments );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void SetCursorToCenterOfVolume ( tkm_tVolumeType iVolume ) {

  xVoxel    anaIdx;
  tkm_tErr  eResult = tkm_tErr_NoErr;
  MWin_tErr eWindow = MWin_tErr_NoErr;

  DebugEnterFunction( ("SetCursorToCenterOfVolume ( iVolume=%d )",
                       (int)iVolume) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* Get the size of the volume and set the cursor. */
  DebugNote( ("Sending color scale update to tcl window") );

  /* Set a voxel to the center of the volume. */
  xVoxl_Set( &anaIdx,
             gnAnatomicalDimensionX/2,
             gnAnatomicalDimensionY/2,
             gnAnatomicalDimensionZ/2 );

  /* Tell the window to go there. */
  eWindow = MWin_SetCursor( gMeditWindow, -1, &anaIdx );
  DebugAssertThrowX( (MWin_tErr_NoErr == eWindow),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* Zoom around new cursor. */
  eWindow = MWin_SetZoomCenterToCursor( gMeditWindow, -1 );
  DebugAssertThrowX( (MWin_tErr_NoErr == eWindow),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* Zoom out. */
  eWindow = MWin_SetZoomLevel( gMeditWindow, -1, 1 );
  DebugAssertThrowX( (MWin_tErr_NoErr == eWindow),
                     eResult, tkm_tErr_ErrorAccessingVolume );

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
                                        xVoxelRef       iaMRIIdx,
                                        int             inCount,
                                        Volm_tValue     inLow,
                                        Volm_tValue     inHigh,
                                        Volm_tValue     inNewValue ) {

  int          nVoxel   = 0;
  float        value    = 0;

  if ( NULL == gAnatomicalVolume[iVolume] ) {
    return;
  }

  /* for each voxel we got... */
  for ( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {

    /* get the value at this point. */
    Volm_GetValueAtMRIIdx( gAnatomicalVolume[iVolume],
                           &(iaMRIIdx[nVoxel]), &value );

    /* if it's in the range and it's a different value... */
    if ( value >= inLow && value <= inHigh &&
         value != inNewValue ) {

      Volm_SetValueAtMRIIdx( gAnatomicalVolume[iVolume],
                             &(iaMRIIdx[nVoxel]), inNewValue );

      switch ( iVolume ) {
      case tkm_tVolumeType_Main:
        /* if this is an edit on the main volume, add it to the undo list
           using the EditAnatomicalVolume function and also add it to the
           undo volume. */
        AddVoxelAndValueToUndoList( EditAnatomicalVolume,
                                    &(iaMRIIdx[nVoxel]), value );
        AddMRIIdxAndValueToUndoVolume( &(iaMRIIdx[nVoxel]), value );
        break;

      case tkm_tVolumeType_Aux:
        /* if this is an edit on the aux volume, add it to the undo list
           using the EditAuxAnatomicalVolume function. */
        AddVoxelAndValueToUndoList( EditAuxAnatomicalVolume,
                                    &(iaMRIIdx[nVoxel]), value );
        break;

      default:
        break;
      }
    }
  }

  /* volume is dirty. */
  SetVolumeDirty( iVolume, TRUE );
}

void CloneAnatomicalVolumeInRangeArray ( tkm_tVolumeType iDestVolume,
                                         tkm_tVolumeType iSourceVolume,
                                         xVoxelRef       iaMRIIdx,
                                         int             inCount,
                                         Volm_tValue     inLow,
                                         Volm_tValue     inHigh ) {

  int          nVoxel   = 0;
  float        sourceValue    = 0;
  float        value    = 0;

  if ( NULL == gAnatomicalVolume[iDestVolume] ||
       NULL == gAnatomicalVolume[iSourceVolume] ) {
    return;
  }

  /* for each voxel we got... */
  for ( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {

    /* get the value at this point. */
    Volm_GetValueAtMRIIdx( gAnatomicalVolume[iDestVolume],
                           &(iaMRIIdx[nVoxel]), &value );

    /* Get the new value from the source volume. */
    Volm_GetValueAtMRIIdx( gAnatomicalVolume[iSourceVolume],
                           &(iaMRIIdx[nVoxel]), &sourceValue );

    /* if it's in the range and it's a different value... */
    if ( value >= inLow && value <= inHigh &&
         value != sourceValue ) {

      Volm_SetValueAtMRIIdx( gAnatomicalVolume[iDestVolume],
                             &(iaMRIIdx[nVoxel]), sourceValue );

      switch ( iDestVolume ) {
      case tkm_tVolumeType_Main:
        /* if this is an edit on the main volume, add it to the undo list
           using the EditAnatomicalVolume function and also add it to the
           undo volume. */
        AddVoxelAndValueToUndoList( EditAnatomicalVolume,
                                    &(iaMRIIdx[nVoxel]), value );
        AddMRIIdxAndValueToUndoVolume( &(iaMRIIdx[nVoxel]), value );
        break;

      case tkm_tVolumeType_Aux:
        /* if this is an edit on the aux volume, add it to the undo list
           using the EditAuxAnatomicalVolume function. */
        AddVoxelAndValueToUndoList( EditAuxAnatomicalVolume,
                                    &(iaMRIIdx[nVoxel]), value );
        break;

      default:
        break;
      }
    }
  }

  /* volume is dirty. */
  SetVolumeDirty( iDestVolume, TRUE );
}

void SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
                                 int             iMRIIdxX0,
                                 int             iMRIIdxX1,
                                 int             iMRIIdxY0,
                                 int             iMRIIdxY1,
                                 int             iMRIIdxZ0,
                                 int             iMRIIdxZ1,
                                 float           iNewValue ) {

  tkm_tErr   eResult = tkm_tErr_NoErr;
  Volm_tErr  eVolume = Volm_tErr_NoErr;
  xVoxel     beginMRIIdx;
  xVoxel     endMRIIdx;
  xVoxel     curMRIIdx;

  DebugEnterFunction( ("SetAnatomicalVolumeRegion( iVolume=%d, "
                       "iMRIIdxX0=%d, iMRIIdxY0=%d, iMRIIdxZ0=%d, "
                       "iMRIIdxX1=%d, iMRIIdxY1=%d, iMRIIdxZ1=%d )",
                       iVolume, iMRIIdxX0, iMRIIdxY0, iMRIIdxZ0,
                       iMRIIdxX1, iMRIIdxY1, iMRIIdxZ1) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume < tkm_knNumVolumeTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* Make sure we got a good range. */
  xVoxl_Set( &beginMRIIdx,
             MIN( iMRIIdxX0, iMRIIdxX1 ),
             MIN( iMRIIdxY0, iMRIIdxY1 ),
             MIN( iMRIIdxZ0, iMRIIdxZ1 ) );
  eVolume =
    Volm_VerifyIdxInMRIBounds( gAnatomicalVolume[iVolume], &beginMRIIdx );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_InvalidParameter );
  xVoxl_Set( &endMRIIdx,
             MAX( iMRIIdxX0, iMRIIdxX1 ),
             MAX( iMRIIdxY0, iMRIIdxY1 ),
             MAX( iMRIIdxZ0, iMRIIdxZ1 ) );
  eVolume =
    Volm_VerifyIdxInMRIBounds( gAnatomicalVolume[iVolume], &endMRIIdx );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_InvalidParameter );

  /* step through the volume and set everything to this value. */
  DebugNote( ("Setting volume values") );
  xVoxl_Copy( &curMRIIdx, &beginMRIIdx );
  while ( xVoxl_IncrementWithMinsUntilLimits( &curMRIIdx,
                                              xVoxl_GetX(&beginMRIIdx),
                                              xVoxl_GetY(&beginMRIIdx),
                                              xVoxl_GetX(&endMRIIdx),
                                              xVoxl_GetY(&endMRIIdx),
                                              xVoxl_GetZ(&endMRIIdx) ) ) {
    Volm_SetValueAtMRIIdx_( gAnatomicalVolume[iVolume], &
                            curMRIIdx, iNewValue );
  }

  /* volume is dirty. */
  SetVolumeDirty( iVolume, TRUE );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

int EditAnatomicalVolume ( xVoxelRef iMRIIdx, int inValue ) {

  tkm_tErr    eResult = tkm_tErr_NoErr;
  Volm_tErr   eVolume = Volm_tErr_NoErr;
  float       value   = 0;
  char        sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("EditAnatomicalVolume( iMRIIdx=%p, inValue=%d)",
                       iMRIIdx, inValue) );

  /* get the current value so we can return it. set the new value */
  eVolume = Volm_GetValueAtMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                   iMRIIdx, &value );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );
  eVolume = Volm_SetValueAtMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Main],
                                   iMRIIdx, inValue );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* resend the min/max values if necessary. */
  if ( inValue < gfaAnaColorMin[tkm_tVolumeType_Main] ||
       inValue > gfaAnaColorMax[tkm_tVolumeType_Main] ) {

    Volm_GetValueMinMax( gAnatomicalVolume[tkm_tVolumeType_Main],
                         &gfaAnaColorMin[tkm_tVolumeType_Main],
                         &gfaAnaColorMax[tkm_tVolumeType_Main] );

    sprintf( sTclArguments, "%d %.2f %.2f", (int)tkm_tVolumeType_Main,
             gfaAnaColorMin[tkm_tVolumeType_Main],
             gfaAnaColorMax[tkm_tVolumeType_Main] );

    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValueMinMax,
                        sTclArguments );
  }


  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return (int)value;
}

int EditAuxAnatomicalVolume ( xVoxelRef iMRIIdx, int inValue ) {

  tkm_tErr    eResult = tkm_tErr_NoErr;
  Volm_tErr   eVolume = Volm_tErr_NoErr;
  float       value   = 0;
  char        sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction(("EditAuxAnatomicalVolume( iMRIIdx=%p, inValue=%d)",
                      iMRIIdx, inValue) );

  /* get the current value so we can return it. set the new value */
  eVolume = Volm_GetValueAtMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Aux],
                                   iMRIIdx, &value );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );
  eVolume = Volm_SetValueAtMRIIdx( gAnatomicalVolume[tkm_tVolumeType_Aux],
                                   iMRIIdx, inValue );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* resend the min/max values if necessary. */
  if ( inValue < gfaAnaColorMin[tkm_tVolumeType_Aux] ||
       inValue > gfaAnaColorMax[tkm_tVolumeType_Aux] ) {

    Volm_GetValueMinMax( gAnatomicalVolume[tkm_tVolumeType_Aux],
                         &gfaAnaColorMin[tkm_tVolumeType_Aux],
                         &gfaAnaColorMax[tkm_tVolumeType_Aux] );

    sprintf( sTclArguments, "%d %.2f %.2f", (int)tkm_tVolumeType_Aux,
             gfaAnaColorMin[tkm_tVolumeType_Aux],
             gfaAnaColorMax[tkm_tVolumeType_Aux] );

    tkm_SendTclCommand( tkm_tTclCommand_UpdateVolumeValueMinMax,
                        sTclArguments );
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return (int)value;
}

tkm_tErr FloodFillAnatomicalVolume ( tkm_tVolumeType iVolume,
                                     xVoxelRef       iAnaIdx,
                                     int             inValue,
                                     tBoolean        ib3D,
                                     float           iFuzzy,
                                     float           iDistance ) {

  tkm_tErr                    eResult      = tkm_tErr_NoErr;
  Volm_tFloodParams           params;
  tkm_tFloodFillAnatomicalCallbackData  callbackData;
  Volm_tErr                   eVolume      = Volm_tErr_NoErr;
  xVoxel                      mriIdx;

  DebugEnterFunction( ("FloodFillAnatomicalVolume( iVolume=%d, "
                       "iAnaIdx=%d,%d,%d, iVnalue=%d, ib3D=%d, iFuzzy=%f "
                       "iDistance=%f )", iVolume, xVoxl_ExpandInt(iAnaIdx),
                       inValue, ib3D, iFuzzy, iDistance) );

  DebugAssertThrowX( (NULL != iAnaIdx), eResult, tkm_tErr_InvalidParameter );

  /* We need an MRI index. */
  Volm_ConvertIdxToMRIIdx( gAnatomicalVolume[iVolume], iAnaIdx, &mriIdx );

  xVoxl_Copy( &params.mSourceIdx, &mriIdx );
  params.mfFuzziness             = iFuzzy;
  params.mComparatorType         = Volm_tValueComparator_EQ;
  params.mComparatorFunc         = NULL;
  params.mfMaxDistance           = iDistance;
  params.mb3D                    = ib3D;
  MWin_GetOrientation ( gMeditWindow, &params.mOrientation );

  /* Set some local parameters telling us what to do on the
     callback. */
  callbackData.mnValue  = inValue;
  callbackData.mnCount  = 0;
  callbackData.mVolume  = iVolume;

  /* Set the callback function data. Tell it to use the callback data
     we just initialized. */
  params.mpFunction     = FloodFillAnatomicalCallback;
  params.mpFunctionData = (void*)&callbackData;

  /* Now get the source volume value. */
  Volm_GetValueAtIdx( gAnatomicalVolume[iVolume],
                      iAnaIdx, &params.mfSourceValue );

  /* Start listening for a cancel. */
  xUtil_StartListeningForUserCancel();

  /* Do it! */
  eVolume = Volm_Flood( gAnatomicalVolume[iVolume], &params );

  /* If we filled more than 1000 voxels, we printed a message and
     started printing update dots. Now close off the message. */
  if ( callbackData.mnCount > 1000 ) {
    printf( "done. %d voxels filled. \n", callbackData.mnCount );
  }

  /* Stop listening for the cancel. */
  xUtil_StopListeningForUserCancel();

  UpdateAndRedraw();

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tVisitCommand FloodFillAnatomicalCallback ( xVoxelRef iMRIIdx,
                                                 float     iValue,
                                                 void*     iData ) {
  tkm_tFloodFillAnatomicalCallbackData* callbackData;

  /* Grab our callback data and edit the volume. */
  callbackData = (tkm_tFloodFillAnatomicalCallbackData*)iData;
  if ( tkm_tVolumeType_Main == callbackData->mVolume ) {

    EditAnatomicalVolume( iMRIIdx, callbackData->mnValue );

    /* if this is an edit on the main volume, add it to the undo list
       using the EditAnatomicalVolume function and also add it to the
       undo volume. */
    AddVoxelAndValueToUndoList( EditAnatomicalVolume, iMRIIdx, iValue );

    AddMRIIdxAndValueToUndoVolume( iMRIIdx, iValue );

  } else {

    EditAuxAnatomicalVolume( iMRIIdx, callbackData->mnValue );

    /* if this is an edit on the aux volume, add it to the undo list
       using the EditAuxAnatomicalVolume function. */
    AddVoxelAndValueToUndoList( EditAuxAnatomicalVolume, iMRIIdx, iValue );
  }

  /* Incremenet our count. If it's over 1000, print a message saying
     the user can cancel and start printing update dots. */
  callbackData->mnCount++;
  if ( callbackData->mnCount == 1000 ) {
    printf( "Filling (press ctrl-c to cancel) " );
  }
  if ( callbackData->mnCount > 1000 &&
       callbackData->mnCount % 100 == 0 ) {
    printf( "." );
    fflush( stdout );
  }

  /* Check the user cancel. If they canceled, stop. */
  if ( xUtil_DidUserCancel() ) {
    return Volm_tVisitComm_Stop;
  }

  return Volm_tVisitComm_Continue;
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

tkm_tErr LoadFunctionalOverlay( char* isFileName,
                                char* isOffsetFileName,
                                FunD_tRegistrationType iRegType,
                                char* isRegistrationFileName ) {

  tkm_tErr  eResult                               = tkm_tErr_NoErr;
  char      sFileName[tkm_knPathLen]              = "";
  char      sOffsetFileName[tkm_knPathLen]        = "";
  char      sRegistrationFileName[tkm_knPathLen]  = "";
  char*     psOffsetFileName                      = NULL;
  char*     psRegistrationFileName                = NULL;
  FunV_tErr eFunctional                           = FunV_tErr_NoError;
  char      sError[tkm_knErrStringLen]            = "";

  DebugEnterFunction( ("LoadFunctionalOverlay( isFileName=%s, "
                       "isOffsetFileName=%s, iRegType=%d, "
                       "isRegistrationFileName=%s )", isFileName,
                       isOffsetFileName, iRegType , isRegistrationFileName) );

  /* Make our filename. If we have an offset or registration file
     name, make that too. */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_PWD,
                sFileName, sizeof(sFileName) );

  if ( NULL != isOffsetFileName ) {
    DebugNote( ("Making file name from %s", isOffsetFileName) );
    MakeFileName( isOffsetFileName, tkm_tFileName_PWD,
                  sOffsetFileName, sizeof(sOffsetFileName) );
    psOffsetFileName = sOffsetFileName;
  } else {
    psOffsetFileName = NULL;
  }

  if ( NULL != isRegistrationFileName ) {
    DebugNote( ("Making file name from %s", isRegistrationFileName) );
    MakeFileName( isRegistrationFileName, tkm_tFileName_PWD,
                  sRegistrationFileName, sizeof(sRegistrationFileName) );
    psRegistrationFileName = sRegistrationFileName;
  } else {
    psRegistrationFileName = NULL;
  }

  /* Load the overlay. */
  DebugNote( ("Loading overlay") );
  eFunctional = FunV_LoadOverlay( gFunctionalVolume,
                                  sFileName, psOffsetFileName,
                                  iRegType, psRegistrationFileName,
                                  gAnatomicalVolume[tkm_tVolumeType_Main]);
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
                     eResult, tkm_tErr_CouldntLoadOverlay );

  /* turn overlay display on */
  MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_FunctionalOverlay,
                       TRUE );

  if ( gMeditWindow ) {
    tkm_SendTclCommand( tkm_tTclCommand_ShowFuncOverlayOptions, "1" );
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  xUtil_snprintf( sError, sizeof(sError), "Loading functional overlay %s",
                  isFileName );
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

tkm_tErr LoadFunctionalTimeCourse( char* isFileName,
                                   char* isOffsetFileName,
                                   FunD_tRegistrationType iRegType,
                                   char* isRegistrationFileName ) {

  tkm_tErr  eResult                               = tkm_tErr_NoErr;
  char      sFileName[tkm_knPathLen]              = "";
  char      sOffsetFileName[tkm_knPathLen]        = "";
  char      sRegistrationFileName[tkm_knPathLen]  = "";
  char*     psOffsetFileName                      = NULL;
  char*     psRegistrationFileName                = NULL;
  FunV_tErr eFunctional                           = FunV_tErr_NoError;
  char      sError[tkm_knErrStringLen]            = "";

  DebugEnterFunction( ("LoadFunctionalTimeCourse( isFileName=%s, "
                       "isOffsetFileName=%s, iRegType=%d, "
                       "isRegistrationFileName=%s )", isFileName,
                       isOffsetFileName, iRegType, isRegistrationFileName) );

  /* Make our filename. If we have an offset or registration file
     name, make that too. */
  DebugNote( ("Making file name from %s", isFileName) );
  MakeFileName( isFileName, tkm_tFileName_PWD,
                sFileName, sizeof(sFileName) );

  if ( NULL != isOffsetFileName ) {
    DebugNote( ("Making file name from %s", isOffsetFileName) );
    MakeFileName( isOffsetFileName, tkm_tFileName_PWD,
                  sOffsetFileName, sizeof(sOffsetFileName) );
    psOffsetFileName = sOffsetFileName;
  } else {
    psOffsetFileName = NULL;
  }

  if ( NULL != isRegistrationFileName ) {
    DebugNote( ("Making file name from %s", isRegistrationFileName) );
    MakeFileName( isRegistrationFileName, tkm_tFileName_PWD,
                  sRegistrationFileName, sizeof(sRegistrationFileName) );
    psRegistrationFileName = sRegistrationFileName;
  } else {
    psRegistrationFileName = NULL;
  }

  /* Load the overlay. */
  DebugNote( ("Loading overlay") );
  eFunctional = FunV_LoadTimeCourse( gFunctionalVolume,
                                     sFileName, psOffsetFileName,
                                     iRegType, psRegistrationFileName,
                                     gAnatomicalVolume[tkm_tVolumeType_Main] );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
                     eResult, tkm_tErr_CouldntLoadOverlay );

  /* turn overlay display on */
  MWin_SetDisplayFlag( gMeditWindow, -1, DspA_tDisplayFlag_FunctionalOverlay,
                       TRUE );


  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  xUtil_snprintf( sError, sizeof(sError), "Loading functional overlay %s",
                  isFileName );
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

  switch ( isAxis ) {
  case 'x':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'y':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'z':
    switch ( orientation ) {
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
    default:
      goto cleanup;
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

  switch ( isDirection ) {
  case 'x':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'y':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'z':
    switch ( orientation ) {
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
    default:
      goto cleanup;
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

  switch ( isDirection ) {
  case 'x':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'y':
    switch ( orientation ) {
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
    default:
      goto cleanup;
      break;
    }
    break;
  case 'z':
    switch ( orientation ) {
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
    default:
      goto cleanup;
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

  DebugEnterFunction( ("NewSegmentationVolume( iVolume=%d, iFromVolume=%d, "
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
  eVolume = Volm_CreateFromVolume( newVolume, gAnatomicalVolume[iFromVolume]);
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  eVolume = Volm_SetAllValues( newVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Scale up if we need to. */
  if ( gbScaleUpVolume ) {
    DebugNote( ("Setting min voxel size on new seg volume to 1") );
    eVolume = Volm_SetMinVoxelSizeToOne( newVolume );
    DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }

  // see if there is a ctab stored in volume first
  if (newVolume->mpMriValues->ct != NULL)
  {
    gColorTable[iVolume] = newVolume->mpMriValues->ct ;
    SendColorTableInformationToTcl( iVolume );
  }
  else
  {
    /* Try to load the color table. */
    DebugNote( ("Loading color table.") );
    eResult = LoadSegmentationColorTable( iVolume, isColorFileName );
    DebugAssertThrowX( (tkm_tErr_NoErr == eResult),
                       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }


  /* allocate flag volume from the new segmentation volume. set
     everything in it to zero */
  DebugNote( ("Creating segmentation flag volume") );
  eVolume = Volm_DeepClone( newVolume, &newChangedVolume );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_CouldntLoadSegmentation );
  eVolume = Volm_SetAllValues( newChangedVolume, 0 );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );



  /* free existing segmentations if present. */
  if ( NULL != gSegmentationVolume[iVolume] ) {
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
  if ( gMeditWindow ) {
    if ( tkm_tSegType_Main == iVolume ) {
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

  if ( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError), "Creating Segmentation" );
    tkm_DisplayError( sError, tkm_GetErrorString(eResult),
                      "Tkmedit couldn't create a segmentation." );
  }

  if ( NULL != newVolume ) {
    Volm_Delete( &newVolume );
  }
  if ( NULL != newChangedVolume ) {
    Volm_Delete( &newChangedVolume );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

tkm_tErr LoadSegmentationVolume ( tkm_tSegType iVolume,
                                  char*        isVolumeDirWithPrefix,
                                  char*        isColorFileName,
                                  tBoolean     ibConform ) {

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

  DebugNote( ("Importing segmentation") );
  eVolume = Volm_ImportData( newVolume, sSegmentationFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_CouldntLoadSegmentation );
  if ( gbScaleUpVolume ) {
    Volm_SetMinVoxelSizeToOne( newVolume );
  }

  /* Conform if desired. */
  if ( ibConform ) {
    Volm_Conform( newVolume );
  }

  if (newVolume->mpMriValues->ct != NULL) // ctab embedded in volume
  {
    gColorTable[iVolume] = newVolume->mpMriValues->ct ;
    SendColorTableInformationToTcl( iVolume );
  }
  else /* Try to load the color table. */
  {
    DebugNote( ("Loading color table.") );
    eResult = LoadSegmentationColorTable( iVolume, isColorFileName );
    DebugAssertThrowX( (Volm_tErr_NoErr == eResult),
                       eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  }

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
  if ( NULL != gSegmentationVolume[iVolume] ) {
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
  if ( gMeditWindow ) {
    if ( tkm_tSegType_Main == iVolume ) {
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

  if ( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError),
                    "Loading Segmentation %s", isVolumeDirWithPrefix );
    tkm_DisplayError( sError,
                      tkm_GetErrorString(eResult),
                      "Tkmedit couldn't read the segmentation you "
                      "specified. This could be because the segmentation "
                      "volume wasn't a valid COR volume directory." );
  }

  if ( NULL != newVolume ) {
    Volm_Delete( &newVolume );
  }
  if ( NULL != newChangedVolume ) {
    Volm_Delete( &newChangedVolume );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

void SaveSegmentationVolume ( tkm_tSegType iVolume,
                              char*        isFileName,
                              tBoolean     ibSetFileName ) {

  tkm_tErr  eResult         = tkm_tErr_NoErr;
  char      sError[tkm_knErrStringLen] = "";
  Volm_tErr eVolume         = Volm_tErr_NoErr;
  char      sSegmentationFileName[tkm_knPathLen] = "";
  char*     psSegmentationFileName     = NULL;

  DebugEnterFunction( ("SaveSegmentationVolume ( iVolume=%d, "
                       "isFileName=%s )", (int)iVolume, isFileName) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
                     eResult, tkm_tErr_InvalidParameter );
  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume]),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* make a file name if they gave us one. */
  if ( NULL != isFileName ) {
    MakeFileName( isFileName, tkm_tFileName_Segmentation,
                  sSegmentationFileName, sizeof(sSegmentationFileName) );
    psSegmentationFileName = sSegmentationFileName;
  } else {
    psSegmentationFileName = NULL;
  }

  /* Export the thing to a COR volume. */
  OutputPrint "Saving... " EndOutputPrint;
  DebugNote( ("Exporting to COR") );
  eVolume = Volm_Save( gSegmentationVolume[iVolume],
                       psSegmentationFileName, ibSetFileName );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_CouldntWriteFile );
  OutputPrint "done.\n" EndOutputPrint;

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  if ( eResult != tkm_tErr_CouldntLoadColorTable ) {
    xUtil_snprintf( sError, sizeof(sError),
                    "Saving Segmentation %s", isFileName );
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
  if ( NULL != isVolumeDir ) {
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
  eVolume = Volm_Save( newVolume, sSegmentationFileName, FALSE );
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
  tBoolean     bInternalTable          = FALSE;
  COLOR_TABLE* colorTable              = NULL;
  int          nVertex                 = 0;
  VERTEX*      pVertex                 = 0;
  int          nRed                    = 0;
  int          nGreen                  = 0;
  int          nBlue                   = 0;
  int          nStructure              = 0;
  xVoxel       surfRAS;
  xVoxel       MRIIdx;
  float        dx                      = 0;
  float        dy                      = 0;
  float        dz                      = 0;
  float        len                     = 0;
  float        d                       = 0;
  char         sError[tkm_knErrStringLen] = "";

  DebugEnterFunction( ("ImportSurfaceAnnotationToSegmentation( iVolume=%d "
                       "isAnnotationFileName=%s, isColorFileName=%s )",
                       (int)iVolume, isAnnotationFileName, isColorFileName) );

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

  /* Annotations usually have embedded color tables. If so, the
     surface will have a link to it now. Check if there's one
     present. If so... */
  Surf_IsInternalColorTablePresent( gSurface[tkm_tSurfaceType_Main],
                                    &bInternalTable );
  if ( bInternalTable ) {

    /* Make a new color table from the embedded one. If we got
       it... */
    eSurface = Surf_NewColorTableFromInternal( gSurface[tkm_tSurfaceType_Main],
                                               &colorTable );
    if ( Surf_tErr_NoErr == eSurface ) {

      /* Replace our existing color table with it. */
      CTABfree( &gColorTable[iVolume] );
      gColorTable[iVolume] = colorTable;

      /* Update the window and tcl stuff. */
      MWin_SetSegmentationColorTable( gMeditWindow, iVolume,
                                      -1, gColorTable[iVolume] );

      SendColorTableInformationToTcl( iVolume );
    }
  }

  /* For every vertex in the surface, get the annot value. Unpack it
     to an RGB value. Look in the LUT to get the corresponding
     index. Get the orig coords from the vertex and set the
     segmentation value at those coords to the index we got. */
  fprintf( stdout, "Converting annotation..." );
  for ( nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
    pVertex = &mris->vertices[nVertex];

    /* Get the color and then the index. */
    if ( 0 != pVertex->annotation ) {
      MRISAnnotToRGB( pVertex->annotation, nRed, nGreen, nBlue );
      CTABfindRGBi( gColorTable[iVolume], nRed, nGreen, nBlue, &nStructure );

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

      for ( d = 0 ; d <= len; d = d+0.1 ) {
        xVoxl_SetFloat( &surfRAS,
                        pVertex->x + (d * dx),
                        pVertex->y + (d * dy),
                        pVertex->z + (d * dz) );
        if ( gbUseRealRAS ) {
          eVolume =
            Volm_ConvertRASToMRIIdx( newVolume, &surfRAS, &MRIIdx );
        } else {
          eVolume =
            Volm_ConvertSurfaceRASToMRIIdx( newVolume, &surfRAS, &MRIIdx );
        }
        eVolume =
          Volm_SetValueAtMRIIdx( newVolume, &MRIIdx, (float)nStructure );
      }
    }
    if ( !(nVertex % 1000) ) {
      fprintf( stdout, "\rConverting annotation... %.2f%% done",
               ((float)nVertex / (float)mris->nvertices) * 100.0 );
      fflush( stdout );
    }
  }
  fprintf( stdout, "\rConverting annotation... 100%% done       \n" );


  /* free existing segmentations if present. */
  if ( NULL != gSegmentationVolume[iVolume] ) {
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
  if ( gMeditWindow ) {
    if ( tkm_tSegType_Main == iVolume ) {
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

  if ( eResult != tkm_tErr_CouldntLoadColorTable ) {

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
  gColorTable[iVolume] = CTABreadASCII( isColorFileName );
  DebugAssertThrowX( (NULL != gColorTable[iVolume]),
                     eResult, tkm_tErr_CouldntLoadColorTable );

  /* Send the labels to tcl. */
  SendColorTableInformationToTcl( iVolume );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  xUtil_snprintf( sError, sizeof(sError),
                  "Loading Color Table %s", isColorFileName );
  tkm_DisplayError( sError,
                    tkm_GetErrorString(eResult),
                    "Tkmedit couldn't read the color table you "
                    "specified. This could be because the file wasn't "
                    "valid or wasn't found." );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

tkm_tErr SendColorTableInformationToTcl ( tkm_tSegType iVolume ) {

  tkm_tErr  eResult                        = tkm_tErr_NoErr;
  int       eCTAB                          = 0;
  int       nNumEntries                    = 0;
  int       nEntry                         = 0;
  int       bValid                         = FALSE;
  char      sLabel[1024]                   = "";
  char      sTclArguments[tkm_knTclCmdLen] = "";
  char      sFileName[256]                 = "";

  DebugEnterFunction( ("SendColorTableInformationToTcl( iVolume=%d )",
                       iVolume) );

  /* build the color table for the interface. go through the entries
     and for each one, make a string of its number and its label, and send
     that as a new entry. */
  DebugNote( ("Clearing color table in interface") );
  tkm_SendTclCommand( tkm_tTclCommand_ClearSegColorTable, "" );

  /* Find out how many total entries we have. */
  DebugNote( ("Getting number of total color table entries") );
  CTABgetNumberOfTotalEntries( gColorTable[iVolume], &nNumEntries );

  /* Iterate over all the entries, but only get names for the valid
     ones. */
  for ( nEntry = 0; nEntry < nNumEntries; nEntry++ ) {

    /* If not valid, skip it. */
    CTABisEntryValid( gColorTable[iVolume], nEntry, &bValid );
    if (!bValid)
      continue;

    DebugNote( ("Getting label for entry %d/%d", nEntry, nNumEntries) );
    eCTAB =
      CTABcopyName( gColorTable[iVolume], nEntry, sLabel, sizeof(sLabel));
    DebugAssertThrowX( (NO_ERROR == eCTAB),
                       eResult, tkm_tErr_Unrecoverable );
    if ( strcmp( sLabel, "" ) == 0 )
      xUtil_strncpy( sLabel, "None", sizeof(sLabel) );

    DebugNote( ("Making tcl command") );
    xUtil_snprintf( sTclArguments, sizeof(sTclArguments),
                    "%d \"%s\"", nEntry, sLabel );
    tkm_SendTclCommand( tkm_tTclCommand_AddSegColorTableEntry,
                        sTclArguments );
  }

  /* Copy the file name and send that. */
  CTABcopyFileName( gColorTable[iVolume], sFileName, sizeof(sFileName) );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateSegmentationColorTable,
                      sFileName );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


void SetSegmentationAlpha ( float ifAlpha ) {

  char sTclArguments[STRLEN] = "";

  if ( ifAlpha >= 0 && ifAlpha <= 1.0 ) {

    gfSegmentationAlpha = ifAlpha;

    if ( NULL != gMeditWindow ) {
      MWin_SetSegmentationAlpha( gMeditWindow, -1, ifAlpha );
      UpdateAndRedraw ();
    }

  }

  sprintf( sTclArguments, "%f", gfSegmentationAlpha );
  tkm_SendTclCommand( tkm_tTclCommand_UpdateSegmentationVolumeAlpha,
                      sTclArguments );

}

void GetSegmentationColorAtVoxel ( tkm_tSegType iVolume,
                                   xVoxelRef    iMRIIdx,
                                   xColor3fRef  iBaseColor,
                                   xColor3fRef  oColor ) {

  int        nIndex       = 0;
  int        eCTAB        = NO_ERROR;
  float      fRed         = 0;
  float      fGreen       = 0;
  float      fBlue        = 0;
  float      fAlpha       = 0;
  float      fFinalAlpha  = 0;

  /* get the index of this voxel */
  GetSegLabel( iVolume, iMRIIdx, &nIndex, NULL );

  /* If 0, just use the base color */
  if ( 0 == nIndex ) {
    *oColor = *iBaseColor;
    goto cleanup;
  }

  /* get the color out of the color map */
  eCTAB = CTABrgbaAtIndexf( gColorTable[iVolume], nIndex,
                            &fRed, &fGreen, &fBlue, &fAlpha );
  if ( NO_ERROR != eCTAB )
    goto error;

  fFinalAlpha = gfSegmentationAlpha * fAlpha;

  /* blend them */
  oColor->mfRed = (fFinalAlpha * fRed) +
    (float)((1.0-fFinalAlpha) * iBaseColor->mfRed);
  oColor->mfGreen = (fFinalAlpha * fGreen) +
    (float)((1.0-fFinalAlpha) * iBaseColor->mfGreen);
  oColor->mfBlue = (fFinalAlpha * fBlue) +
    (float)((1.0-fFinalAlpha) * iBaseColor->mfBlue);

  goto cleanup;

 error:

  DebugPrint( ( "Error in GetSegmentationColor.\n" ) );

  *oColor = *iBaseColor;

 cleanup:
  return;
}

void GetSegLabel ( tkm_tSegType iVolume,
                   xVoxelRef    iMRIIdx,
                   int*         onIndex,
                   char*        osLabel ) {

  int         eCTAB        = NO_ERROR;
  Volm_tErr   eVolume      = Volm_tErr_NoErr;
  int         index        = 0;
  float       fValue       = 0;

  /* get the voxel at this location */
  eVolume =
    Volm_GetValueAtMRIIdx( gSegmentationVolume[iVolume], iMRIIdx, &fValue );
  if ( Volm_tErr_NoErr == eVolume ) {
    index = (int)fValue;
  } else {
    index = 0;
  }

  /* get the label out of the table and return it and in the index */
  if ( NULL != osLabel ) {
    if ( 0 == index ) {
      strcpy( osLabel, "None" );
    } else {
      eCTAB = CTABcopyName( gColorTable[iVolume], index, osLabel, 256 );
      if ( NO_ERROR != eCTAB ) {
        /* pass an out of bounds notice. */
        strcpy( osLabel, "Out of bounds." );
      }
    }
  }

  /* return the index */
  if ( NULL != onIndex ) {
    *onIndex = index;
  }

  goto cleanup;

  goto error;
 error:

  DebugPrint( ( "Error in GetSegmentationColor.\n" ) );

  if ( NULL != onIndex )
    *onIndex = -1;
  if ( NULL != osLabel )
    strcpy( osLabel, "None" );

 cleanup:
  return;
}

Volm_tVisitCommand AddSimilarVoxelToSelection ( xVoxelRef iMRIIdx,
                                                float     iValue,
                                                void*     ipnTarget ) {
  int nIndex = 0;
  int nTargetIndex = 0;

  nIndex = (int) iValue;
  nTargetIndex = *(int*)ipnTarget;

  if ( nIndex == nTargetIndex )
    AddVoxelsToSelection( iMRIIdx, 1 );

  return Volm_tVisitComm_Continue;
}

Volm_tVisitCommand AddSimilarVoxelToGraphAvg ( xVoxelRef iMRIIdx,
                                               float     iValue,
                                               void*     ipnTarget ) {
  int    nIndex       = 0;
  int    nTargetIndex = 0;

  nIndex = (int) iValue;
  nTargetIndex = *(int*)ipnTarget;

  if ( nIndex == nTargetIndex ) {
    FunV_AddMRIIdxToSelectionRange( gFunctionalVolume, iMRIIdx );
  }

  return Volm_tVisitComm_Continue;
}

Volm_tVisitCommand SetChangedSegmentationValue ( xVoxelRef iMRIIdx,
                                                 float     iValue,
                                                 void*     iData ) {

  tkm_tExportSegmentationParamsRef params = NULL;
  float                            value;

  if ( iValue ) {

    params = (tkm_tExportSegmentationParamsRef)iData;

    /* Get the value from the segmentation volume and set it in the
       volume to export. */
    Volm_GetValueAtIdx( params->mSrcVolume, iMRIIdx, &value );
    Volm_SetValueAtIdx( params->mDestVolume, iMRIIdx, value );
  }

  return Volm_tVisitComm_Continue;
}


tkm_tErr SelectSegLabelAtCursor ( tkm_tSegType iVolume ) {

  tkm_tErr  eResult     = tkm_tErr_NoErr;
  Volm_tErr eVolume     = Volm_tErr_NoErr;
  xVoxel    cursorIdx;
  float     value       = 0;

  DebugEnterFunction( ("SelectSegLabel ( iVolume=%d )", iVolume) );

  DebugAssertThrowX( (iVolume >= 0 && iVolume <= tkm_knNumSegTypes),
                     eResult, tkm_tErr_InvalidParameter );

  /* make sure we're loaded */
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume] ||
                      NULL != gColorTable[iVolume]),
                     eResult, tkm_tErr_SegmentationNotLoaded );

  /* Get the cursor */
  MWin_GetCursor( gMeditWindow, &cursorIdx );

  /* Get the value at the cursor. */
  eVolume = Volm_GetValueAtIdx( gSegmentationVolume[iVolume],
                                &cursorIdx, &value );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  /* Select this value. */
  SelectSegLabel( iVolume, value );

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;

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
  DebugAssertThrowX( (NULL != gSegmentationVolume[iVolume] ||
                      NULL != gColorTable[iVolume]),
                     eResult, tkm_tErr_SegmentationNotLoaded );

  /* check entry index */
  DebugNote( ("Getting number of entries") );
  CTABgetNumberOfTotalEntries( gColorTable[iVolume], &nNumEntries );
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
  FunV_tErr eFunctional    = FunV_tErr_NoError;
  int       nNumEntries    = 0;
  char      sLabel[1024]   = "";
  sLabel[0]=0;

  DebugEnterFunction( ("GraphSegLabel ( iVolume=%d, inIndex=%d )",
                       iVolume, inIndex) );

  /* make sure we're loaded */
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != gSegmentationVolume ||
                      NULL != gColorTable[iVolume]),
                     eResult, tkm_tErr_SegmentationNotLoaded );

  /* check entry index */
  DebugNote( ("Getting number of entries") );
  CTABgetNumberOfTotalEntries( gColorTable[iVolume], &nNumEntries );
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
                         &AddSimilarVoxelToGraphAvg, (void*)&inIndex );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );
  OutputPrint "done.\n" EndOutputPrint;

  /* finish the selection */
  eFunctional = FunV_EndSelectionRange( gFunctionalVolume );
  DebugAssertThrowX( (FunV_tErr_NoError == eFunctional),
                     eResult, tkm_tErr_ErrorAccessingFunctionalVolume );

  /* set the graph window */
  CTABcopyName( gColorTable[iVolume], inIndex, sLabel, sizeof(sLabel) );
  if ( sLabel[0] != 0 )
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
  if ( gbAcceptingTclCommands && gDisplayIntermediateResults ) {
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
  if ( NULL != gPreviousSegmentationVolume[iVolume] ) {
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
  GCAreclassifyUsingGibbsPriors
    (
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
                       xVoxelRef    iMRIIdx,
                       int          inIndex ) {

  float       fOldValue = 0;
  static int  nVolumeIndexBugs = 0;

  DebugEnterFunction( ("EditSegmentation( iVolume=%d, iMRIIdx=%p, "
                       "inIndex=%d )", iVolume, iMRIIdx, inIndex) );

  /* For some reason, David was getting random bugs where iVolume
     would be way out of bounds. I'm trying to track this down while
     still letting him work. */
  if ( iVolume != tkm_tSegType_Main ) {
    nVolumeIndexBugs++;
    if ( nVolumeIndexBugs == 1 ) {
      fprintf( stderr,
               "ATTENTION: Please send the .xdebug_tkmedit file\n"
               "           to freesurfer@nmr.mgh.harvard.edu\n"
               "           when you're done.\n");
    }
    if ( nVolumeIndexBugs < 5 ) {
      xDbg_Printf( "EditSegmentation: iVolume was %d. Stack:\n", iVolume );
      xDbg_PrintStack ();
    }
    iVolume = tkm_tSegType_Main;
  }

  /* Get the old value for the undo entry. */
  DebugNote( ("Getting old value") );
  Volm_GetValueAtMRIIdx( gSegmentationVolume[iVolume], iMRIIdx, &fOldValue );

  /* Change the value in the segmentation and changed volumes. */
  DebugNote( ("Setting value in segmentation volume") );
  Volm_SetValueAtMRIIdx( gSegmentationVolume[iVolume],
                         iMRIIdx, (float)inIndex );
  DebugNote( ("Setting value in segmentation changed volume") );
  Volm_SetValueAtMRIIdx( gSegmentationChangedVolume[iVolume],
                         iMRIIdx, (float)1 );

  DebugExitFunction;

  return (int)fOldValue;
}

void CalcSegLabelVolume ( tkm_tSegType iVolume,
                          xVoxelRef    iMRIIdx,
                          int*         oCount ) {

  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  Volm_tFloodParams           params;
  int       index   = 0;
  tkm_tSumSimilarSegmentationParams callbackData;

  DebugEnterFunction( ("CalcSegLabelVolume( iMRIIdx=%p, oCount=%p )",
                       iMRIIdx, oCount) );
  DebugAssertThrowX( (NULL != iMRIIdx && NULL != oCount),
                     eResult, tkm_tErr_InvalidParameter );


  DebugNote( ("Getting initial value at %d,%d,%d",
              xVoxl_ExpandInt( iMRIIdx )) );
  GetSegLabel( iVolume, iMRIIdx, &index, NULL );
  DebugAssertQuietThrow( (0 != index) );

  xVoxl_Copy( &params.mSourceIdx, iMRIIdx );
  params.mfSourceValue           = index;
  params.mfFuzziness             = 0;
  params.mComparatorType         = Volm_tValueComparator_EQ;
  params.mComparatorFunc         = NULL;
  params.mfMaxDistance           = pow(256,3);
  params.mb3D                    = TRUE;
  MWin_GetOrientation ( gMeditWindow, &params.mOrientation );

  /* Initialize callback data */
  callbackData.mCount = 0;
  callbackData.mSourceLabel = index;

  /* SEt the callback data. */
  params.mpFunction = SumSimilarValues;
  params.mpFunctionData = (void*)&callbackData;

  /* Run the flood. */
  DebugNote( ("Running compare") );
  eVolume =
    Volm_Flood( gSegmentationVolume[iVolume], &params );
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingSegmentationVolume );

  DebugNote( ("Setting return value") );
  *oCount = callbackData.mCount;

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

Volm_tVisitCommand SumSimilarValues ( xVoxelRef iMRIIdx,
                                      float     iValue,
                                      void*     iData ) {

  tkm_tSumSimilarSegmentationParams*  params = NULL;

  params = (tkm_tSumSimilarSegmentationParams*)iData;

  /* If the values are the same, increment the count. */
  if ( FEQUAL( iValue, params->mSourceLabel ) ) {
    params->mCount++;
  }

  return Volm_tVisitComm_Continue;
}


void SetSegmentationValue ( tkm_tSegType iVolume,
                            xVoxelRef    iMRIIdx,
                            int          inIndex ) {

  DebugEnterFunction( ("SetSegmentationValue( iVolume=%d, iMRIIdx=%p "
                       "inIndex=%d )", iVolume, iMRIIdx, inIndex) );

  DebugNote( ("Passing to EditSegmentation") );
  EditSegmentation( iVolume, iMRIIdx, inIndex );

  DebugExitFunction;
}

void SetSegmentationValues ( tkm_tSegType iVolume,
                             xVoxelRef    iaMRIIdx,
                             int          inCount,
                             int          inIndex ) {

  int nVoxel = 0;

  DebugEnterFunction( ("SetSegmentationValues( iVolume=%d, iaMRIIdx=%p "
                       "inCount=%d, inIndex=%d )", iVolume, iaMRIIdx,
                       inCount, inIndex) );

  for ( nVoxel = 0; nVoxel < inCount; nVoxel++ ) {
    EditSegmentation( iVolume, &(iaMRIIdx[nVoxel]), inIndex );
  }

  DebugExitFunction;
}


tkm_tErr FloodFillSegmentation ( tkm_tSegType      iVolume,
                                 xVoxelRef         iAnaIdx,
                                 int               inIndex,
                                 tBoolean          ib3D,
                                 tkm_tVolumeTarget iSrc,
                                 float             iFuzzy,
                                 float             iDistance ) {

  tkm_tErr                    eResult      = tkm_tErr_NoErr;
  Volm_tFloodParams           params;
  tkm_tFloodFillCallbackData  callbackData;
  mriVolumeRef                sourceVolume = NULL;
  Volm_tErr                   eVolume      = Volm_tErr_NoErr;
  char                        sTclArguments[tkm_knTclCmdLen] = "";

  DebugEnterFunction( ("FloodFillSegmentation( iVolume=%d, iAnaIdx=%d,%d,%d "
                       "inIndex=%d, ib3D=%d, iSrc=%d, iFuzzy=%f "
                       "iDistance=%f", iVolume, xVoxl_ExpandInt(iAnaIdx),
                       inIndex, ib3D, iSrc, iFuzzy, iDistance) );

  DebugAssertThrowX( (NULL != iAnaIdx), eResult, tkm_tErr_InvalidParameter );

  xVoxl_Copy( &params.mSourceIdx, iAnaIdx );
  params.mfFuzziness             = iFuzzy;
  params.mComparatorType         = Volm_tValueComparator_EQ;
  params.mComparatorFunc         = NULL;
  params.mfMaxDistance           = iDistance;
  params.mb3D                    = ib3D;
  MWin_GetOrientation ( gMeditWindow, &params.mOrientation );

  /* Set some local parameters telling us what to do on the
     callback. */
  callbackData.mTargetVolume = iVolume;
  callbackData.mnNewSegLabel = inIndex;
  callbackData.mnCount       = 0;

  /* Set the callback function data. Tell it to use the callback data
     we just initialized. */
  params.mpFunction     = FloodFillSegmentationCallback;
  params.mpFunctionData = (void*)&callbackData;

  /* See what source volume to use. For everything other than the main
     anatomical volume, check to see if it exists first. For the
     segmentation targets, set the fuzzinees to 0, since it doesn't
     really apply to label values. */
  switch ( iSrc ) {
  case tkm_tVolumeTarget_MainAna:
    sourceVolume = gAnatomicalVolume[tkm_tVolumeType_Main];
    break;
  case tkm_tVolumeTarget_AuxAna:
    if ( NULL == gAnatomicalVolume[tkm_tVolumeType_Aux] ) {
      strcpy( sTclArguments,
              "\"Cannot use aux volume as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gAnatomicalVolume[tkm_tVolumeType_Aux];
    break;
  case tkm_tVolumeTarget_MainSeg:
    if ( NULL == gSegmentationVolume[tkm_tSegType_Main] ) {
      strcpy( sTclArguments,
              "\"Cannot use segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gSegmentationVolume[tkm_tSegType_Main];
    params.mfFuzziness = 0;
    break;
  case tkm_tVolumeTarget_AuxSeg:
    if ( NULL == gSegmentationVolume[tkm_tSegType_Aux] ) {
      strcpy
        ( sTclArguments,
          "\"Cannot use aux segmentation as source if it is not loaded.\"" );
      tkm_SendTclCommand( tkm_tTclCommand_ErrorDlog, sTclArguments );
      goto cleanup;
    }
    sourceVolume = gSegmentationVolume[tkm_tSegType_Aux];
    params.mfFuzziness = 0;
    break;
  default:
    DebugThrowX( eResult, tkm_tErr_InvalidParameter );
    break;
  }

  /* Now get the source volume value. */
  Volm_GetValueAtIdx( sourceVolume, iAnaIdx, &params.mfSourceValue );

  /* Start listening for a cancel. */
  xUtil_StartListeningForUserCancel();

  /* Do it! */
  eVolume = Volm_Flood( sourceVolume, &params );

  /* If we selected more than 1000 voxels, we printed a message and
     started printing update dots. Now close off the message. */
  if ( callbackData.mnCount > 1000 ) {
    printf( "done. %d voxels filled. \n", callbackData.mnCount );
  }

  /* Stop listening for the cancel. */
  xUtil_StopListeningForUserCancel();

  UpdateAndRedraw();

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tVisitCommand FloodFillSegmentationCallback ( xVoxelRef iAnaIdx,
                                                   float     iValue,
                                                   void*     iData ) {
  tkm_tFloodFillCallbackData* callbackData;

  /* Grab our callback data and edit the segmentation. */
  callbackData = (tkm_tFloodFillCallbackData*)iData;
  EditSegmentation( callbackData->mTargetVolume, iAnaIdx,
                    callbackData->mnNewSegLabel );

  /* Incremenet our count. If it's over 1000, print a message saying
     the user can cancel and start printing update dots. */
  callbackData->mnCount++;
  if ( callbackData->mnCount == 1000 ) {
    printf( "Filling (press ctrl-c to cancel) " );
  }
  if ( callbackData->mnCount > 1000 &&
       callbackData->mnCount % 100 == 0 ) {
    printf( "." );
    fflush( stdout );
  }

  /* Check the user cancel. If they canceled, stop. */
  if ( xUtil_DidUserCancel() ) {
    return Volm_tVisitComm_Stop;
  }

  return Volm_tVisitComm_Continue;
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
  if ( Volm_tErr_NoErr != eVolm ) {

    /* Try it again with a BFLOAT hint. */
    DebugNote( ("Copying BFLOAT hint to string") );
    sprintf( sNameWithHint, "%s@BFLOAT", isNameEV );

    DebugNote( ("Importing data %s", sNameWithHint) );
    eVolm = Volm_ImportData( EVVolume, sNameWithHint );
    if ( Volm_tErr_NoErr != eVolm ) {

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
  if ( Volm_tErr_NoErr != eVolm ) {

    DebugNote( ("Copying BFLOAT hint to string") );
    sprintf( sNameWithHint, "%s@BFLOAT", isNameFA );

    DebugNote( ("Importing data %s", sNameWithHint) );
    eVolm = Volm_ImportData( FAVolume, sNameWithHint );
    if ( Volm_tErr_NoErr != eVolm ) {

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

    for ( nFrame = 0; nFrame < 3; nFrame++ ) {
      Volm_GetValueAtIdxFrameUnsafe( EVVolume, &screenIdx, nFrame, &EVValue );
      Volm_SetValueAtIdxFrame( EVVolume, &screenIdx, nFrame,
                               EVValue * MIN( 1, FAValue) );
    }

    if ( xVoxl_GetY( &EVIdx ) == 0 ) {
      fprintf( stdout, "\rProcessing DTI volumes: %.2f%% done",
               (xVoxl_GetFloatZ( &EVIdx ) / (float)zEVZ) * 100.0 );
    }

  } while ( xVoxl_IncrementUntilLimits( &EVIdx, zEVX, zEVY, zEVZ ) );

  fprintf( stdout, "\rProcessing DTI volumes: 100%% done.         \n" );

  /* Delete the FA volume. */
  Volm_Delete( &FAVolume );

  /* If we already have DTI volumes, delete it. */
  if ( NULL != gDTIVolume ) {
    Volm_Delete( &gDTIVolume );
  }

  /* Use these DTI volumes. */
  gDTIVolume = EVVolume;

  /* Set it in the window. */
  if ( NULL != gMeditWindow ) {
    DebugNote( ("Setting DTI volume in main window") );
    MWin_SetDTIVolume( gMeditWindow, -1, gDTIVolume );
  }

  /* Save the axis -> component settings */
  gaDTIAxisForComponent[xColr_tComponent_Red] = iRedAxis;
  gaDTIAxisForComponent[xColr_tComponent_Green] = iGreenAxis;
  gaDTIAxisForComponent[xColr_tComponent_Blue] = iBlueAxis;

  if ( NULL != gMeditWindow ) {
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

void SetDTIAlpha ( float ifAlpha ) {

  char sTclArguments[STRLEN] = "";

  if ( ifAlpha >= 0 && ifAlpha <= 1.0 ) {

    gfDTIAlpha = ifAlpha;

    if ( NULL != gMeditWindow ) {
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
  if ( xUndL_tErr_NoErr != eList ) {
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
  if ( NULL == this ) {
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
  if ( NULL == ioEntry ) {
    DebugPrint( ( "DeleteUndoEntry(): NULL parameter.\n" ) );
    return;
  }

  this = *ioEntry;
  if ( NULL == this ) {
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

  if ( NULL == this ) {
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
  if ( xUndL_tErr_NoErr != eList ) {
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
  if ( NULL == entry ) {
    DebugPrint( ( "AddVoxelAndValueToUndoList(): Couldn't create entry.\n" ) );
    return;
  }

  /* add the entry */
  eList = xUndL_AddEntry( gUndoList, entry );
  if ( xUndL_tErr_NoErr != eList ) {
    DebugPrint( ( "AddVoxelAndValueToUndoList(): Error in xUndL_AddEntry "
                  "%d: %s\n", eList, xUndL_GetErrorString( eList ) ) );
  }
}

void RestoreUndoList () {

  xUndL_tErr eList = xUndL_tErr_NoErr;

  /* restore the list */
  eList = xUndL_Restore( gUndoList );
  if ( xUndL_tErr_NoErr != eList ) {
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

  if ( NULL == *ipEntryToDelete ) {
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
  if ( NULL == entryToUndo ) {
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

tkm_tErr InitUndoVolume () {

  tkm_tErr   eResult = tkm_tErr_NoErr;
  Volm_tErr  eVolume = Volm_tErr_NoErr;
  int        nZ      = 0;
  int        nY      = 0;

  DebugEnterFunction( ("InitUndoVolume()") );

  /* delete the volume if it exists */
  if ( NULL != gUndoVoxelValueVolume ||
       NULL != gUndoVoxelFlagVolume )
    DeleteUndoVolume();

  /* Get dimensions of the main MRI. */
  eVolume = Volm_GetDimensions( gAnatomicalVolume[tkm_tVolumeType_Main],
                                &gUndoVolumeDimensions[0],
                                &gUndoVolumeDimensions[1],
                                &gUndoVolumeDimensions[2]);
  DebugAssertThrowX( (Volm_tErr_NoErr == eVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* allocate the flag volume */
  DebugNote( ("Creating gUndoVoxelFlagVolume") );
  gUndoVoxelFlagVolume =
    (tBoolean***) calloc( gUndoVolumeDimensions[2], sizeof(tBoolean**) );
  DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume),
                     eResult, tkm_tErr_CouldntAllocate );
  for ( nZ = 0; nZ < gUndoVolumeDimensions[2]; nZ++ ) {
    gUndoVoxelFlagVolume[nZ] =
      (tBoolean**) calloc( gUndoVolumeDimensions[1], sizeof(tBoolean*) );
    DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume[nZ]),
                       eResult, tkm_tErr_CouldntAllocate );
    for ( nY = 0; nY < gUndoVolumeDimensions[1]; nY++ ) {
      gUndoVoxelFlagVolume[nZ][nY] =
        (tBoolean*) calloc( gUndoVolumeDimensions[0], sizeof(tBoolean) );
      DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume[nZ][nY]),
                         eResult, tkm_tErr_CouldntAllocate );
    }
  }

  /* allocate the value volume */
  DebugNote( ("Creating gUndoVoxelValueVolume") );
  gUndoVoxelValueVolume =
    (Volm_tValue***) calloc( gUndoVolumeDimensions[2], sizeof(Volm_tValue**) );
  DebugAssertThrowX( (NULL != gUndoVoxelValueVolume),
                     eResult, tkm_tErr_CouldntAllocate );
  for ( nZ = 0; nZ < gUndoVolumeDimensions[2]; nZ++ ) {
    gUndoVoxelValueVolume[nZ] =
      (Volm_tValue**) calloc( gUndoVolumeDimensions[1], sizeof(Volm_tValue*) );
    DebugAssertThrowX( (NULL != gUndoVoxelValueVolume[nZ]),
                       eResult, tkm_tErr_CouldntAllocate );
    for ( nY = 0; nY < gUndoVolumeDimensions[1]; nY++ ) {
      gUndoVoxelValueVolume[nZ][nY] =
        (Volm_tValue*) calloc( gUndoVolumeDimensions[0], sizeof(Volm_tValue) );
      DebugAssertThrowX( (NULL != gUndoVoxelValueVolume[nZ][nY]),
                         eResult, tkm_tErr_CouldntAllocate );
    }
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

void DeleteUndoVolume () {

  int       nZ      = 0;
  int       nY      = 0;

  DebugEnterFunction( ("DeleteUndoVolume()") );

  if ( NULL != gUndoVoxelValueVolume ) {
    DebugNote( ("Deleting gUndoVoxelValueVolume") );
    for ( nZ = 0; nZ < gUndoVolumeDimensions[2]; nZ++ ) {
      for ( nY = 0; nY < gUndoVolumeDimensions[1]; nY++ ) {
        free( gUndoVoxelValueVolume[nZ][nY] );
      }
      free( gUndoVoxelValueVolume[nZ] );
    }
    free( gUndoVoxelValueVolume );
    gUndoVoxelValueVolume = NULL;
  }

  if ( NULL != gUndoVoxelFlagVolume ) {
    DebugNote( ("Deleting gUndoVoxelFlagVolume") );
    for ( nZ = 0; nZ < gUndoVolumeDimensions[2]; nZ++ ) {
      for ( nY = 0; nY < gUndoVolumeDimensions[1]; nY++ ) {
        free( gUndoVoxelFlagVolume[nZ][nY] );
      }
      free( gUndoVoxelFlagVolume[nZ] );
    }
    free( gUndoVoxelFlagVolume );
    gUndoVoxelFlagVolume = NULL;
  }

  DebugExitFunction;
}

void AddMRIIdxAndValueToUndoVolume ( xVoxelRef    iMRIIdx,
                                     Volm_tValue  iValue ) {

  tkm_tErr eResult = tkm_tErr_NoErr;

  DebugEnterFunction( ("AddMRIIdxAndValueToUndoVolume( iMRIIdx=%p, iValue=%f)",
                       iMRIIdx, iValue) );

  DebugNote( ("Checking if undo volumes exist") );
  DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume &&
                      NULL != gUndoVoxelValueVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  DebugNote( ("Checking params") );
  DebugAssertThrowX( (NULL != iMRIIdx),
                     eResult, tkm_tErr_InvalidParameter );

  /* Set the flag here. */
  gUndoVoxelFlagVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)] = 1;

  /* Set the value here. */
  gUndoVoxelValueVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)] = iValue;

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

tBoolean IsMRIIdxInUndoVolume ( xVoxelRef iMRIIdx ) {

  tkm_tErr eResult     = tkm_tErr_NoErr;
  tBoolean bIsInVolume = FALSE;

  DebugEnterFunction( ("IsMRIIdxInUndoVolume( iMRIIdx=%p )", iMRIIdx) );

  DebugNote( ("Checking if undo volumes exist") );
  DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume &&
                      NULL != gUndoVoxelValueVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  DebugNote( ("Checking params") );
  DebugAssertThrowX( (NULL != iMRIIdx),
                     eResult, tkm_tErr_InvalidParameter );

  /* Return the flag here. */
  bIsInVolume =
    gUndoVoxelFlagVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)];

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return bIsInVolume;
}

void RestoreUndoVolumeAroundMRIIdx ( xVoxelRef iMRIIdx ) {

  tkm_tErr  eResult = tkm_tErr_NoErr;
  Volm_tErr eVolume = Volm_tErr_NoErr;
  int       nZ      = 0;
  int       nY      = 0;
  int       nX      = 0;
  xVoxel    MRIIdx;


  DebugEnterFunction( ("IsMRIIdxInUndoVolume( iMRIIdx=%p )", iMRIIdx) );

  DebugNote( ("Checking if undo volumes exist") );
  DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume &&
                      NULL != gUndoVoxelValueVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  DebugNote( ("Checking params") );
  DebugAssertThrowX( (NULL != iMRIIdx),
                     eResult, tkm_tErr_InvalidParameter );


  /* if this voxel is in the volume... */
  if ( gUndoVoxelFlagVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)] ) {

    /* Restore the value */
    Volm_SetValueAtMRIIdx_
      ( gAnatomicalVolume[tkm_tVolumeType_Main], iMRIIdx,
        (Volm_tValue)gUndoVoxelValueVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)] );

    /* "Remove" it from the volume. */
    gUndoVoxelFlagVolume[xVoxl_GetZ(iMRIIdx)][xVoxl_GetY(iMRIIdx)][xVoxl_GetX(iMRIIdx)] = FALSE;

    /* Try restoring surrounding voxels as well. */
    for ( nZ = xVoxl_GetZ(iMRIIdx)-1; nZ <= xVoxl_GetZ(iMRIIdx)+1; nZ++ )
      for ( nY = xVoxl_GetY(iMRIIdx)-1; nY <= xVoxl_GetY(iMRIIdx)+1; nY++ )
        for ( nX = xVoxl_GetX(iMRIIdx)-1; nX <= xVoxl_GetX(iMRIIdx)+1; nX++ ) {
          xVoxl_Set( &MRIIdx, nX, nY, nZ );
          eVolume =
            Volm_VerifyMRIIdx_( gAnatomicalVolume[tkm_tVolumeType_Main],
                                &MRIIdx);
          if ( Volm_tErr_NoErr == eVolume ) {
            RestoreUndoVolumeAroundMRIIdx( &MRIIdx );
          }
        }

    UpdateAndRedraw();
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
}

void ClearUndoVolume () {

  tkm_tErr  eResult = tkm_tErr_NoErr;
  int       nX      = 0;
  int       nY      = 0;
  int       nZ      = 0;

  DebugEnterFunction( ("DeleteUndoVolume()") );

  DebugNote( ("Checking if undo volumes exist") );
  DebugAssertThrowX( (NULL != gUndoVoxelFlagVolume &&
                      NULL != gUndoVoxelValueVolume),
                     eResult, tkm_tErr_ErrorAccessingVolume );

  /* Set all the flags to FALSE. */
  DebugNote( ("Setting flags in gUndoVoxelFlagVolume to false") );
  for ( nZ = 0; nZ < gUndoVolumeDimensions[2]; nZ++ ) {
    for ( nY = 0; nY < gUndoVolumeDimensions[1]; nY++ ) {
      for ( nX = 0; nX < gUndoVolumeDimensions[0]; nX++ ) {
        gUndoVoxelFlagVolume[nZ][nY][nX] = FALSE;
      }
    }
  }

  DebugCatch;
  DebugCatchError( eResult, tkm_tErr_NoErr, tkm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;
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
  MATRIX*   anaIdxToRASTransform = NULL;
  MATRIX*   RAStoAnaIdxTransform = NULL;

  DebugEnterFunction( ("LoadHeadPts( isHeadPtsFile=%s, isTransformFile=%s )",
                       isHeadPtsFile, isTransformFile) );

  /* make filenames */
  DebugNote( ("Making file name from %s", isHeadPtsFile) );
  MakeFileName( isHeadPtsFile, tkm_tFileName_HeadPoints,
                sHeadPtsFile, sizeof(sHeadPtsFile) );
  if ( NULL != isTransformFile ) {
    DebugNote( ("Making file name from %s", isTransformFile) );
    MakeFileName( isTransformFile, tkm_tFileName_HeadPoints,
                  sTransformFile, sizeof(sTransformFile) );
    spTransformFileArg = sTransformFile;
  } else {
    spTransformFileArg = NULL;
  }

  /* We want surface RAS -> ana Idx for our client transform, as the
     head points are in surface RAS space. */
  DebugNote( ("Getting surface RAS->anaIdx transform from volume") );
  RAStoAnaIdxTransform =
    voxelFromSurfaceRAS_(gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues);
  DebugAssertThrowX( (NULL != RAStoAnaIdxTransform),
                     eResult, tkm_tErr_ErrorAccessingTransform );

  /* Read the head points. */
  DebugNote( ("Creating head points list") );
  eHeadPts = HPtL_New( &gHeadPoints, sHeadPtsFile, spTransformFileArg,
                       RAStoAnaIdxTransform );
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

  if( NULL != anaIdxToRASTransform ) {
    DebugNote( ("Freeing anaIdxToRASTransform") );
    MatrixFree( &anaIdxToRASTransform );
  }

  if( NULL != RAStoAnaIdxTransform ) {
    DebugNote( ("Freeing RAStoAnaIdxTransform") );
    MatrixFree( &RAStoAnaIdxTransform );
  }

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

  switch ( isAxis ) {
  case 'x':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_Y );
      break;
    default:
      goto cleanup;
      break;
    }
    break;
  case 'y':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    default:
      goto cleanup;
      break;
    }
    break;
  case 'z':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_Y );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Rotate( gHeadPoints, ifDegrees, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Rotate( gHeadPoints, -ifDegrees, tAxis_X );
      break;
    default:
      goto cleanup;
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

  switch ( isDirection ) {
  case 'x':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_X );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Y );
      break;
    default:
      goto cleanup;
      break;
    }
    break;
  case 'y':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints, ifDistance, tAxis_Z );
      break;
    default:
      goto cleanup;
      break;
    }
    break;
  case 'z':
    switch ( orientation ) {
    case mri_tOrientation_Coronal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Y );
      break;
    case mri_tOrientation_Horizontal:
      HPtL_Translate( gHeadPoints, -ifDistance, tAxis_Z );
      break;
    case mri_tOrientation_Sagittal:
      HPtL_Translate( gHeadPoints,  -ifDistance, tAxis_X );
      break;
    default:
      goto cleanup;
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
  if ( NULL == isNewLabel )
    goto error;

  /* get the point from the window */
  eWindow = MWin_GetSelectedHeadPt( gMeditWindow, &pHeadPoint );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  /* if we got one... */
  if ( NULL != pHeadPoint ) {

    /* set the name */
    strcpy( pHeadPoint->msLabel, isNewLabel );
  }

  goto cleanup;

 error:

  if ( MWin_tErr_NoErr != eWindow ) {
    DebugPrint( ( "MWin error %d in SetSelectedHeadPointLabel: %s\n",
                  eWindow, MWin_GetErrorString( eWindow ) ) );
  }

 cleanup:

  return;
}

void AlignSelectedHeadPointToMRIIdx ( xVoxelRef iMRIIdx ) {

  HPtL_tHeadPointRef pHeadPoint = NULL;
  MWin_tErr       eWindow  = MWin_tErr_NoErr;
  HPtL_tErr       eHeadPts  = HPtL_tErr_NoErr;

  /* get the point from the window */
  eWindow = MWin_GetSelectedHeadPt( gMeditWindow, &pHeadPoint );
  if ( MWin_tErr_NoErr != eWindow )
    goto error;

  /* if we got one... */
  if ( NULL != pHeadPoint ) {

    /* align to it */
    eHeadPts = HPtL_AlignPointToClientVoxel( gHeadPoints,
                                             pHeadPoint, iMRIIdx );
    if ( HPtL_tErr_NoErr != eHeadPts )
      goto error;

    /* redraw */
    UpdateAndRedraw();
  }

  goto cleanup;

 error:

  if ( MWin_tErr_NoErr != eWindow ) {
    DebugPrint( ( "MWin error %d in AlignSelectedHeadPointToMRIIdx: %s\n",
                  eWindow, MWin_GetErrorString( eWindow ) ) );
  }

  if (  HPtL_tErr_NoErr != eHeadPts ) {
    DebugPrint( ( "HPtL error %d in AlignSelectedHeadPointToMRIIdx: %s\n",
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
  if (trans->type == LINEAR_RAS_TO_RAS)
  {
    MRI *mri_tmp = GCAbuildMostLikelyVolume(gca, NULL) ;
    TransformRas2Vox(trans,
                     gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues,
                     mri_tmp) ;
    MRIfree(&mri_tmp) ;
  }

  TransformInvert
    (trans, gAnatomicalVolume[tkm_tVolumeType_Main]->mpMriValues) ;

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

  if ( NULL != gca )
    GCAfree( &gca );
#if 0
  GCAhistoScaleImageIntensities
    (gGCAVolume,
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

  if ( NULL != gGCAVolume )
    GCAfree( &gGCAVolume );
  if ( NULL != gGCATransform )
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

  if ( NULL != gVLI1 )
    VLfree( &gVLI1 );
  if ( NULL != gVLI2 )
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
  if ( NULL != psFileName ) {
    xUtil_strncpy( gsTkTimerFileName, psFileName, sizeof(gsTkTimerFileName) );
  }

  /* open file and count how many valid entries in it */
  nNumEntries = 0;
  DebugNote( ("Opening file %s for counting", gsTkTimerFileName) );
  file = fopen( gsTkTimerFileName, "r" );
  if ( NULL != file ) {
    while ( !feof( file ) ) {
      DebugNote( ("Looking for string and number in file") );
      eSystem = fscanf( file, "%*s %*d" );
      if ( EOF != eSystem )
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
  for ( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Allocating name string %d", nCurEntry) );
    asNames[nCurEntry] = (char*) malloc( sizeof(char) * tkm_knNameLen );
    DebugAssertThrowX( (NULL != asNames[nCurEntry]),
                       eResult, tkm_tErr_CouldntAllocate );
  }

  /* read in all the entries from the file */
  nCurEntry = 0;
  DebugNote( ("Opening file %s for reading", gsTkTimerFileName) );
  file = fopen( gsTkTimerFileName, "r" );
  if ( NULL != file ) {
    while ( !feof( file ) ) {
      DebugNote( ("Getting line %d", nCurEntry) );
      eSystem = fscanf( file, "%s %d",
                        asNames[nCurEntry], &(anTimes[nCurEntry]) );
      if ( 2 == eSystem ) {
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
  for ( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Looking for %s in entry %d", sSubjectName, nCurEntry) );
    if ( strcmp( asNames[nCurEntry], sSubjectName ) == 0 ) {
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
  for ( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
    DebugNote( ("Writing entry %d to file", nCurEntry) );
    fprintf( file, "%s %d\n", asNames[nCurEntry], anTimes[nCurEntry] );
  }
  if ( !bFound ) {
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

  if ( NULL != asNames ) {
    for ( nCurEntry = 0; nCurEntry < nNumEntries; nCurEntry++ ) {
      DebugNote( ("Freeing name %d", nCurEntry) );
      free( asNames[nCurEntry] );
    }
    DebugNote( ("Freeing name array") );
    free( asNames );
  }
  if ( NULL != anTimes ) {
    DebugNote( ("Freeing times array") );
    free( anTimes );
  }
  if ( NULL != file ) {
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

  fprintf( stdout, "%s", isMessage );
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

void tkm_MakeFileName ( char*         isInput,
                        tkm_tFileName iType,
                        char*         osCompleteFileName,
                        int           inDestSize ) {

  MakeFileName( isInput, iType, osCompleteFileName, inDestSize );
}

void tkm_MakeControlPoint ( xVoxelRef iMRIIdx ) {

  if ( NULL == iMRIIdx ) {
    DebugPrint( ( "tkm_MakeControlPoint(): Passed NULL voxel.\n" ) );
    return;
  }

  AddMRIIdxControlPoint( iMRIIdx, gbSavingControlPoints );
}

void tkm_RemoveControlPointWithinDist ( xVoxelRef        iMRIIdx,
                                        mri_tOrientation iPlane,
                                        int              inDistance ) {

  float         fDistance = 0;
  xVoxelRef     pCtrlPt   = NULL;

  /* find the closest control point */
  fDistance = FindNearestMRIIdxControlPoint( iMRIIdx, iPlane, &pCtrlPt );

  /* if we found one and it's in range... */
  if ( NULL != pCtrlPt &&
       fDistance >= 0.0 &&
       fDistance <= (float)inDistance ) {

    /* delete it */
    DeleteMRIIdxControlPoint( pCtrlPt, gbSavingControlPoints );
  }
}

void tkm_WriteControlFile () {

  WriteControlPointFile( );
}

void tkm_EditAnatomicalVolumeInRange( tkm_tVolumeType  iVolume,
                                      xVoxelRef        inVolumeVox,
                                      Volm_tValue      inLow,
                                      Volm_tValue      inHigh,
                                      Volm_tValue      inNewValue ) {

  EditAnatomicalVolumeInRangeArray( iVolume, inVolumeVox, 1,
                                    inLow, inHigh, inNewValue );

}

void tkm_EditAnatomicalVolumeInRangeArray( tkm_tVolumeType  iVolume,
                                           xVoxelRef        iaVolumeVox,
                                           int              inCount,
                                           Volm_tValue      inLow,
                                           Volm_tValue      inHigh,
                                           Volm_tValue      inNewValue ) {

  EditAnatomicalVolumeInRangeArray( iVolume, iaVolumeVox, inCount,
                                    inLow, inHigh, inNewValue );

}

void tkm_CloneAnatomicalVolumeInRangeArray( tkm_tVolumeType  iDestVolume,
                                            tkm_tVolumeType  iSourceVolume,
                                            xVoxelRef        iaVolumeVox,
                                            int              inCount,
                                            Volm_tValue      inLow,
                                            Volm_tValue      inHigh ) {

  CloneAnatomicalVolumeInRangeArray( iDestVolume, iSourceVolume,
                                     iaVolumeVox, inCount,
                                     inLow, inHigh );

}

void tkm_FloodFillAnatomicalVolume ( tkm_tSegType    iVolume,
                                     xVoxelRef       iAnaIdx,
                                     int             inIndex,
                                     tBoolean        ib3D,
                                     float           iFuzzy,
                                     float           iDistance ) {

  FloodFillAnatomicalVolume( iVolume, iAnaIdx, inIndex, ib3D,
                             iFuzzy, iDistance );
}

void tkm_SetVolumeBrightnessContrast ( tkm_tVolumeType iVolume,
                                       float ifBrightness, float ifContrast ) {

  SetVolumeBrightnessAndContrast( iVolume, ifBrightness, ifContrast );
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

tBoolean tkm_IsSelectionPresent () {

  return (gSelectionCount > 0);
}

void tkm_FloodSelect ( xVoxelRef         iSeedAnaIdx,
                       tBoolean          ib3D,
                       tkm_tVolumeTarget iSrc,
                       float             iFuzzy,
                       float             iDistance,
                       tBoolean          ibSelect ) {

  FloodSelect( iSeedAnaIdx, ib3D, iSrc, iFuzzy, iDistance, ibSelect );
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

void tkm_RestoreUndoVolumeAroundMRIIdx ( xVoxelRef iMRIIdx ) {

  RestoreUndoVolumeAroundMRIIdx( iMRIIdx );
}

tBoolean tkm_IsMRIIdxInUndoVolume ( xVoxelRef iMRIIdx ) {

  return IsMRIIdxInUndoVolume( iMRIIdx );
}

void tkm_GetHeadPoint ( xVoxelRef           iMRIIdx,
                        mri_tOrientation    iOrientation,
                        tBoolean            ibFlat,
                        HPtL_tHeadPointRef* opPoint ) {

  HPtL_tErr eHeadPts = HPtL_tErr_NoErr;
  HPtL_tHeadPointRef pPoint = NULL;
  HPtL_tIterationPlane plane = HPtL_tIterationPlane_X;

  if ( gHeadPoints ) {

    switch ( iOrientation ) {
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

    if ( ibFlat ) {
      eHeadPts = HPtL_FindFlattenedNearestPoint( gHeadPoints, plane,
                                                 iMRIIdx, &pPoint );
    } else {
      eHeadPts = HPtL_FindNearestPoint( gHeadPoints, plane,
                                        1.0, iMRIIdx, &pPoint );
    }
    if ( HPtL_tErr_NoErr != eHeadPts )
      goto error;

    if ( NULL == pPoint )
      goto error;

    *opPoint = pPoint;

  }

  goto cleanup;

 error:

  *opPoint = NULL;

 cleanup:
  return;
}

void tkm_WriteVoxelToEditFile ( xVoxelRef inVolumeVox ) {

  WriteVoxelToEditFile ( inVolumeVox );
}

void tkm_ReadCursorFromEditFile () {

  GotoSavedCursor();
}

void tkm_GetAnatomicalVolume ( tkm_tVolumeType iVolume,
                               mriVolumeRef*   opVolume ) {

  if ( iVolume >= 0 && iVolume < tkm_knNumVolumeTypes ) {
    *opVolume = gAnatomicalVolume[iVolume];
  } else {
    opVolume = NULL;
  }
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

void tkm_CalcSegLabelVolume ( tkm_tSegType iVolume,
                              xVoxelRef    iMRIIdx,
                              int*         onVolume ) {

  CalcSegLabelVolume( iVolume, iMRIIdx, onVolume );
}

void tkm_EditSegmentation ( tkm_tSegType iVolume,
                            xVoxelRef    iMRIIdx,
                            int          inIndex ) {

  DebugEnterFunction( ("tkm_EditSegmentation( iVolume=%d, iMRIIdx=%p, "
                       "inIndex=%d )", iVolume, iMRIIdx, inIndex) );

  SetSegmentationValue( iVolume, iMRIIdx, inIndex );

  DebugExitFunction;
}

void tkm_EditSegmentationArray ( tkm_tSegType iVolume,
                                 xVoxelRef    iaMRIIdx,
                                 int          inCount,
                                 int          inIndex ) {

  DebugEnterFunction( ("tkm_EditSegmentationArray( iVolume=%d, iaMRIIdx=%p "
                       "inCount=%d, inIndex=%d )", iVolume, iaMRIIdx,
                       inCount, inIndex) );

  SetSegmentationValues( iVolume, iaMRIIdx, inCount, inIndex );

  DebugExitFunction;
}


void tkm_FloodFillSegmentation ( tkm_tSegType      iVolume,
                                 xVoxelRef         iAnaIdx,
                                 int               inIndex,
                                 tBoolean          ib3D,
                                 tkm_tVolumeTarget iSrc,
                                 float             iFuzzy,
                                 float             iDistance ) {

  FloodFillSegmentation( iVolume, iAnaIdx, inIndex, ib3D, iSrc,
                         iFuzzy, iDistance );
}

void tkm_SetSurfaceDistance    ( xVoxelRef iAnaIdx,
                                 float     ifDistance ) {

  if ( NULL == gSurface[tkm_tSurfaceType_Main] ) {
    return;
  }

  /* This is right to be using ana idx instead of MRI idx because the
     client space for the surface is ana idx (screen space). */
  Surf_SetVertexValue( gSurface[tkm_tSurfaceType_Main],
                       Surf_tVertexSet_Main, Surf_tValueSet_Val,
                       iAnaIdx, ifDistance );
}

void tkm_SetMRIValueInSurface ( xVoxelRef        iAnaIdx,
                                Surf_tVertexSet  iVertexSet,
                                float            ifValue ) {

  if ( NULL == gSurface[tkm_tSurfaceType_Main] ) {
    return;
  }

  /* This is right to be using ana idx instead of MRI idx because the
     client space for the surface is ana idx (screen space). */
  Surf_SetVertexValue( gSurface[tkm_tSurfaceType_Main],
                       iVertexSet, Surf_tValueSet_Val,
                       iAnaIdx, ifValue );
}

void tkm_ShowNearestSurfaceVertex ( Surf_tVertexSet iVertexSet ) {

  if ( iVertexSet < 0 || iVertexSet >= Surf_knNumVertexSets ) {
    return;
  }

  FindNearestSurfaceVertex( iVertexSet );
}

tBoolean tkm_UseRealRAS () {

  return gbUseRealRAS;
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
  "UpdateLineLength",
  "UpdateZoomLevel",
  "UpdateOrientation",
  "UpdateDisplayFlag",
  "UpdateTool",
  "UpdateBrushTarget",
  "UpdateBrushShape",
  "UpdateBrushInfo",
  "UpdateAnatomicalFillInfo",
  "UpdateFloodSelectParams",
  "UpdateCursorColor",
  "UpdateCursorShape",
  "UpdateSurfaceLineWidth",
  "UpdateSurfaceLineColor",
  "UpdateUseRealRAS",
  "UpdateSegBrushInfo",
  "UpdateVolumeColorScaleInfo",
  "UpdateSegmentationVolumeAlpha",
  "UpdateDTIVolumeAlpha",
  "UpdateTimerStatus",
  "UpdateSubjectDirectory",
  "UpdateSegmentationColorTable",
  "UpdateVolumeDirty",
  "UpdateAuxVolumeDirty",
  "UpdateVolumeValueMinMax",
  "UpdateVolumeSampleType",
  "UpdateVolumeResampleMethod",
  "UpdateSurfaceHemi",
  "UpdateVolumeIsConformed",

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
  "tkm_SetEnableGroupStatus tMenuGroup_AuxVolumeOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_DirtyAnatomicalVolume",
  "tkm_SetEnableGroupStatus tMenuGroup_DirtyAuxAnatomicalVolume",
  "tkm_SetEnableGroupStatus tMenuGroup_VolumeMainTransformLoadedOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_VolumeAuxTransformLoadedOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_OverlayOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_TimeCourseOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_SurfaceLoading",
  "tkm_SetEnableGroupStatus tMenuGroup_SurfaceViewing",
  "tkm_SetEnableGroupStatus tMenuGroup_OriginalSurfaceViewing",
  "tkm_SetEnableGroupStatus tMenuGroup_PialSurfaceViewing",
  "tkm_SetEnableGroupStatus tMenuGroup_HeadPoints",
  "tkm_SetEnableGroupStatus tMenuGroup_VLIOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_GCAOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_DTIOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_Registration",
  "tkm_SetEnableGroupStatus tMenuGroup_SegmentationOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_AuxSegmentationOptions",
  "tkm_SetEnableGroupStatus tMenuGroup_GDFLoaded",
  "tkm_SetEnableGroupStatus tMenuGroup_ControlPoints",
  "ClearSegColorTable",
  "AddSegColorTableEntry",

  /* histogram */
  "BarChart_Draw",

  /* interface configuration */
  "wm geometry . ",
  "wm deiconify .; raise .",
  "CsurfInterface",
  "tkm_Finish",

  /* misc */
  "DoResolveUseRealRASDlog",
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
  char sCommand[tkm_knTclCmdLen] = "";
  char* pTitle = NULL;
  char* pWhile = NULL;
  char* pMessage = NULL;
  char* pTmp = NULL;
  char sTclCmd[tkm_knTclCmdLen];
  xGArr_tErr eList = xGArr_tErr_NoErr;


  /* send all our cached commands */
  /* send all our cached commands */
  eList = xGArr_ResetIterator( gCachedTclCommands );
  if ( xGArr_tErr_NoErr != eList )
    goto error;
  while ( (eList = xGArr_NextItem( gCachedTclCommands, (void*)&sTclCmd ))
          == xGArr_tErr_NoErr ) {

    /* if this is an error command */
    if ( strstr( sTclCmd,
                 kTclCommands[tkm_tTclCommand_ErrorDlog] ) != NULL ) {

      /* get the command, go to the first space, copy that into a
         string, and print it to the shell. we parse out the parts of
         the error message command string to make it look pretty. */
      strcpy( sCommand, sTclCmd );

      pTitle = strchr( sCommand, (int)'\"' );
      if ( NULL != pTitle ) {

        pTitle++;

        pTmp = strchr( pTitle, (int)'\"' );
        if ( NULL != pTmp ) {

          *pTmp = '\0';
          pTmp++;

          pWhile = strchr( pTmp, (int)'\"' );
          if ( NULL != pWhile ) {

            pWhile++;

            pTmp = strchr( pWhile, (int)'\"' );
            if ( NULL != pTmp ) {

              *pTmp = '\0';
              pTmp++;

              pMessage = strchr( pTmp, (int)'\"' );
              if ( NULL != pMessage ) {

                pMessage++;

                pTmp = strchr( pMessage, (int)'\"' );
                if ( NULL != pTmp ) {

                  *pTmp = '\0';
                  pTmp++;

                  printf( "\n\n  Error: %s\n", pTitle );
                  printf( "\n  %s\n", pWhile );
                  printf( "\n  %s\n\n", pMessage );
                  fflush( stdout );
                }
              }
            }
          }
        }
      }
    }
  }

  goto cleanup;
 error:

  /* print error message */
  if ( xGArr_tErr_NoErr != eList ) {
    DebugPrint( ("Error %d in PrintCachedTclErrorDlogsToShell: %s\n",
                 eList, xGArr_GetErrorString(eList) ) );
  }
 cleanup:

  return;
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

  for ( nScript = 0; nScript < gnNumCachedScripts; nScript++ ) {

    DebugAssertThrowX( (gsCachedScriptNames[nScript] != NULL),
                       eResult, tkm_tErr_InvalidScriptName );

    DebugNote( ("Parsing script file %s", gsCachedScriptNames[nScript]) );
    eTcl = Tcl_EvalFile( interp, gsCachedScriptNames[nScript] );
    if ( eTcl != TCL_OK ) {
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
  if ( bDirty ) {
    printf( "    : Attempting to save main volume to /tmp directory...\n" );
    SaveVolume( tkm_tVolumeType_Main, "/tmp", FALSE );
    printf( "    : Data was saved to /tmp.\n" );
    printf( "    :\n" );
  }

  printf( "  : Please send the contents of the file .xdebug_tkmedit\n" );
  printf( "  : that should be in this directory to "
          "freesurfer@nmr.mgh.harvard.edu\n");
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

  if ( xVoxl_IsEqualInt( voxelA, voxelB ) ) {
    eResult = xList_tCompare_Match;
  } else {
    eResult = xList_tCompare_GreaterThan;
  }

  return eResult;
}

tkm_tErr
LoadGCARenormalization(char *renormalization_fname) {
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
  while (cp) {
    nlines++ ;
    cp = fgetl(line, 199, fp) ;
  }
  rewind(fp) ;
  printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
  labels = (int *)calloc(nlines, sizeof(int)) ;
  intensities = (float *)calloc(nlines, sizeof(float)) ;
  cp = fgetl(line, 199, fp) ;
  for (i = 0 ; i < nlines ; i++) {
    sscanf(cp, "%e  %e", &f1, &f2) ;
    labels[i] = (int)f1 ;
    intensities[i] = f2 ;
    cp = fgetl(line, 199, fp) ;
  }
  GCArenormalizeIntensities(gGCAVolume, labels, intensities, nlines) ;
  free(labels) ;
  free(intensities) ;
  return(tkm_tErr_NoErr) ;
}
