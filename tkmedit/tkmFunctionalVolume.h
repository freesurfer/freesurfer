/**
 * @file  tkmFunctionalVolume.h
 * @brief Manages functional volume overlay and timecourse plotting
 *
 * Provides an interface for functional volumes as an overlay or
 * timecourse plot. The same volume can be both, or you can have
 * different volumes for each purpose. Handles loading via the
 * mriFunctionaDataAccess code. Keeps track of the current time point
 * and condition, and sends the Tcl commands to interface with
 * tkm_functional.tcl to draw the graph.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/01 01:41:22 $
 *    $Revision: 1.29 $
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


#ifndef tkmFunctionalVolume_h
#define tkmFunctionalVolume_h

#include "xDebug.h"
#include "xTypes.h"
#include "xVoxel.h"
#include "mriFunctionalDataAccess.h"
#include "xList.h"
#include "tkmedit.h"
#include "tcl.h"

typedef enum {
  FunV_tErr_NoError = 0,
  FunV_tErr_AllocationFailed,
  FunV_tErr_DeletionFailed,
  FunV_tErr_ErrorParsingPathAndStem,
  FunV_tErr_ErrorLoadingVolume,
  FunV_tErr_ErrorAccessingInternalVolume,
  FunV_tErr_ErrorAccessingSelectionList,
  FunV_tErr_ErrorConvertingSecondToTimePoint,
  FunV_tErr_ErrorAllocatingOverlayCache,
  FunV_tErr_ErrorOpeningFile,
  FunV_tErr_OverlayNotLoaded,
  FunV_tErr_TimeCourseNotLoaded,
  FunV_tErr_GraphWindowAlreadyInited,
  FunV_tErr_GraphWindowNotInited,
  FunV_tErr_ErrorFindingScriptTclFile,
  FunV_tErr_ErrorParsingScriptTclFile,
  FunV_tErr_InvalidPointer,
  FunV_tErr_InvalidSignature,
  FunV_tErr_InvalidParameter,
  FunV_tErr_InvalidTimePoint,
  FunV_tErr_InvalidCondition,
  FunV_tErr_InvalidMRIIdx,
  FunV_tErr_InvalidThreshold,
  FunV_tErr_InvalidDisplayFlag,
  FunV_tErr_WrongNumberArgs,
  FunV_tErr_ErrorAccessingAnatomicalVolume,
  FunV_tErr_CouldntLoadBrainMask,
  FunV_tErr_InvalidErrorCode,
  FunV_tErr_knNumErrorCodes
} FunV_tErr;

typedef enum {
  FunV_tDisplayFlag_None = -1,
  FunV_tDisplayFlag_Ol_TruncateNegative = 0,
  FunV_tDisplayFlag_Ol_TruncatePositive,
  FunV_tDisplayFlag_Ol_ReversePhase,
  FunV_tDisplayFlag_Ol_OffsetValues,
  FunV_tDisplayFlag_Ol_IgnoreThreshold,
  FunV_tDisplayFlag_Ol_Grayscale,
  FunV_tDisplayFlag_Ol_Opaque,
  FunV_tDisplayFlag_TC_GraphWindowOpen,
  FunV_tDisplayFlag_TC_OffsetValues,
  FunV_tDisplayFlag_TC_PreStimOffset,
  FunV_knNumDisplayFlags
} FunV_tDisplayFlag;

#define FunV_knFirstOverlayDisplayFlag    FunV_tDisplayFlag_Ol_TruncateNegative
#define FunV_knLastOverlayDisplayFlag     FunV_tDisplayFlag_Ol_IgnoreThreshold
#define FunV_knFirstTimeCourseDisplayFlag FunV_tDisplayFlag_TC_GraphWindowOpen
#define FunV_knLastTimeCourseDisplayFlag  FunV_tDisplayFlag_TC_PreStimOffset

typedef enum {
  FunV_tTclCommand_Ol_DoConfigDlog = 0,
  FunV_tTclCommand_Ol_UpdateNumTimePoints,
  FunV_tTclCommand_Ol_UpdateNumConditions,
  FunV_tTclCommand_Ol_UpdateDisplayFlag,
  FunV_tTclCommand_Ol_UpdateDataName,
  FunV_tTclCommand_Ol_UpdateTimePoint,
  FunV_tTclCommand_Ol_UpdateCondition,
  FunV_tTclCommand_Ol_UpdateThreshold,
  FunV_tTclCommand_Ol_UpdateRange,
  FunV_tTclCommand_Ol_ShowOffsetOptions,
  FunV_tTclCommand_TC_DoConfigDlog,
  FunV_tTclCommand_TC_BeginDrawingGraph,
  FunV_tTclCommand_TC_EndDrawingGraph,
  FunV_tTclCommand_TC_DrawGraph,
  FunV_tTclCommand_TC_ClearGraph,
  FunV_tTclCommand_TC_UpdateNumConditions,
  FunV_tTclCommand_TC_UpdateNumTimePoints,
  FunV_tTclCommand_TC_UpdateTimePoint,
  FunV_tTclCommand_TC_UpdateNumPreStimPoints,
  FunV_tTclCommand_TC_UpdateTimeResolution,
  FunV_tTclCommand_TC_UpdateDisplayFlag,
  FunV_tTclCommand_TC_UpdateDataName,
  FunV_tTclCommand_TC_UpdateLocationName,
  FunV_tTclCommand_TC_UpdateGraphData,
  FunV_tTclCommand_TC_UpdateErrorData,
  FunV_tTclCommand_TC_ShowOffsetOptions,
  FunV_knNumTclCommands
} FunV_tTclCommand;

/* Methods of finding the registration. */
typedef enum {
  FunV_tRegistration_None = -1,
  FunV_tRegistration_File = 0,
  FunV_tRegistration_Find,
  FunV_tRegistration_Identity,
  FunV_knNumRegistrationTypes
} FunV_tRegistrationType;

typedef float FunV_tFunctionalValue;
typedef float* FunV_tOverlayCache;

#define FunV_knColorCacheSize 1000
#define FunV_knMaxFindStatsSteps    10000
#define FunV_knMaxFindStatsDistance 20

/* main class structure */
struct tkmFunctionalVolume {

  tSignature mSignature;

  /* default conversion method */
  FunD_tConversionMethod mDefaultConvMethod;

  /* our volumes.mriFunctionalDataRef is actually tkmFunctionalDataAccess */
  mriFunctionalDataRef  mpOverlayVolume;
  mriFunctionalDataRef  mpTimeCourseVolume;
  mriFunctionalDataRef  mpOverlayOffsetVolume;
  mriFunctionalDataRef  mpTimeCourseOffsetVolume;

  /* overlay cache. this looks just like an anatomical volume. */
  FunV_tOverlayCache mOverlayCache;
  int                manCacheDimensions[3];
  int                manCacheOffsets[3];
  tBoolean           mbUseOverlayCache;
  int                mnCachedTimePoint;
  int                mnCachedCondition;

  /* color cache. */
  float                  mafColorCache[FunV_knColorCacheSize][3];
  FunV_tFunctionalValue  mCachedThresholdMin;
  FunV_tFunctionalValue  mCachedThresholdMid;
  FunV_tFunctionalValue  mCachedThresholdSlope;

  /* current drawing state, what section is being displayed */
  xListRef               mpSelectedVoxels;
  int                    mnTimePoint;
  int                    mnCondition;
  FunV_tFunctionalValue  mThresholdMin;
  FunV_tFunctionalValue  mThresholdMid;
  FunV_tFunctionalValue  mThresholdSlope;
  tBoolean               mabDisplayFlags[FunV_knNumDisplayFlags];
  tBoolean               mbGraphInited;
  tBoolean               mbRegistrationEnabled;

  /* functions to access the outside world */
  void       (*mpOverlayChangedFunction)(void);
  void       (*mpSendTkmeditTclCmdFunction)(tkm_tTclCommand,char*);
  char*      (*mpSendTclCommandFunction)(char*);

  int nRawPlot;
  char *RawPlotFiles[10];
  MRI *RawPlotVols[10];

};
typedef struct tkmFunctionalVolume tkmFunctionalVolume;
typedef tkmFunctionalVolume *tkmFunctionalVolumeRef;

typedef enum {
  FunV_tFindStatsComp_Invalid = -1,
  FunV_tFindStatsComp_GTEoLTE = 0,
  FunV_tFindStatsComp_EQ,
  FunV_tFindStatsComp_GTEThresholdMin,
  FunV_knNumFindStatsComp
} FunV_tFindStatsComp;

#define FunV_kSignature 0x00666000

/* allocator an destructor */
FunV_tErr FunV_New    ( tkmFunctionalVolumeRef* oppVolume,
                        void(*ipOverlayChangedFunction)(void),
                        void(*ipSendTkmeditCmdFunction)(tkm_tTclCommand,char*),
                        char*(*ipSendTclCommandFunction)(char*) );
FunV_tErr FunV_Delete ( tkmFunctionalVolumeRef* ioppVolume );

/* loads the data volumes. The transform should be a->b index->RAS
   (only a->RAS transform will be used). */
FunV_tErr FunV_LoadOverlay    ( tkmFunctionalVolumeRef this,
                                char*                  isPathAndStem,
                                char*                  isOffsetPath,
                                FunV_tRegistrationType iRegistrationType,
                                char*                  isRegistration,
                                mriVolumeRef           iAnatomicalVolume );
FunV_tErr FunV_LoadTimeCourse ( tkmFunctionalVolumeRef this,
                                char*                  isPathAndStem,
                                char*                  isOffsetPath,
                                FunV_tRegistrationType iRegistrationType,
                                char*                  isRegistration,
                                mriVolumeRef           iAnatomicalVolume );

/* this loads a specific volume. will delete one if it already exists. */
FunV_tErr FunV_LoadFunctionalVolume_ ( tkmFunctionalVolumeRef this,
                                       mriFunctionalDataRef*  ioppVolume,
                                       char*                  isFileName,
                                       char*                  isHeaderStem,
                                       FunV_tRegistrationType iRegistrationType,
                                       char*                  isRegPath,
                                       mriVolumeRef          iAnatomicalVolume,
                                       tBoolean               ibPrintErrors );

/* sets conversion method in all volumes */
FunV_tErr FunV_SetConversionMethod ( tkmFunctionalVolumeRef this,
                                     FunD_tConversionMethod iMethod );

/* this takes a functional volume and converts it into an anatomical
   space volume, so it can be indexed with anatomical coords. */
FunV_tErr FunV_InitOverlayCache_ ( tkmFunctionalVolumeRef this );
FunV_tErr FunV_SetOverlayCacheValue_ ( tkmFunctionalVolumeRef this,
                                       xVoxelRef              ipMRIIdx,
                                       FunV_tFunctionalValue  iValue );
FunV_tErr FunV_GetOverlayCacheValue_ ( tkmFunctionalVolumeRef this,
                                       xVoxelRef              ipMRIIdx,
                                       FunV_tFunctionalValue* oValue );
FunV_tErr FunV_UseOverlayCache ( tkmFunctionalVolumeRef this,
                                 tBoolean               ibUseCache );

/* builds a cache table for calculating color values for functional
   values based on the current threshold data */
FunV_tErr FunV_InitColorCache_ ( tkmFunctionalVolumeRef this );

/* get status of loadedness */
FunV_tErr FunV_IsOverlayPresent     ( tkmFunctionalVolumeRef this,
                                      tBoolean*              obIsLoaded );
FunV_tErr FunV_IsTimeCoursePresent  ( tkmFunctionalVolumeRef this,
                                      tBoolean*              obIsLoaded );
FunV_tErr FunV_IsOverlayCacheLoaded ( tkmFunctionalVolumeRef this,
                                      tBoolean*              obIsLoaded );


/* smooths the overlay data */
FunV_tErr FunV_SmoothOverlayData ( tkmFunctionalVolumeRef this,
                                   float                  ifSigma );

/* settors. these check values and if valid, sets internal vars. generates
   proper update msgs for tcl */
FunV_tErr FunV_SetTimeResolution   ( tkmFunctionalVolumeRef this,
                                     float                  inTimeResolution );
FunV_tErr FunV_SetNumPreStimPoints ( tkmFunctionalVolumeRef this,
                                     int                    inNumPoints );
FunV_tErr FunV_SetTimeSecond       ( tkmFunctionalVolumeRef this,
                                     int                    inSecond );
FunV_tErr FunV_SetTimePoint        ( tkmFunctionalVolumeRef this,
                                     int                    inTimePoint );
FunV_tErr FunV_SetCondition        ( tkmFunctionalVolumeRef this,
                                     int                    inCondition );
FunV_tErr FunV_SetThreshold        ( tkmFunctionalVolumeRef this,
                                     FunV_tFunctionalValue  iMin,
                                     FunV_tFunctionalValue  iMid,
                                     FunV_tFunctionalValue  iSlope );
FunV_tErr FunV_SetThresholdUsingFDR( tkmFunctionalVolumeRef this,
                                     float                  iRate,
                                     tBoolean               ibMaskToBrain );
FunV_tErr FunV_SetDisplayFlag      ( tkmFunctionalVolumeRef this,
                                     FunV_tDisplayFlag      iFlag,
                                     tBoolean               iNewValue );
FunV_tErr FunV_EnableRegistration  ( tkmFunctionalVolumeRef this,
                                     tBoolean               iNewValue );

/* moving time point */
FunV_tErr FunV_ChangeTimePointBy   ( tkmFunctionalVolumeRef this,
                                     int                    inDelta );

/* allows functional volume to respond to a click. */
FunV_tErr FunV_MRIIdxClicked ( tkmFunctionalVolumeRef this,
                               xVoxelRef              iMRIIdx );

/* modify the overlay registration */
FunV_tErr FunV_ApplyTransformToOverlay      ( tkmFunctionalVolumeRef this,
    MATRIX*             iTransform );
FunV_tErr FunV_TranslateOverlayRegistration ( tkmFunctionalVolumeRef this,
    float                ifDistance,
    tAxis                  iAxis );
FunV_tErr FunV_RotateOverlayRegistration    ( tkmFunctionalVolumeRef this,
    float                  ifDegrees,
    tAxis                  iAxis,
    xVoxelRef        iCenterMRIIdx );
FunV_tErr FunV_ScaleOverlayRegistration     ( tkmFunctionalVolumeRef this,
    float                  ifFactor,
    tAxis                  iAxis );

/* overlay access */

/* basic accessors to values, based on current plane position if
   applicable */
FunV_tErr FunV_GetValueAtMRIIdx ( tkmFunctionalVolumeRef this,
                                  xVoxelRef               ipVoxel,
                                  tBoolean               iSampled,
                                  FunV_tFunctionalValue* opValue );

/* calculate the rgb values for a color */
FunV_tErr FunV_GetColorForValue ( tkmFunctionalVolumeRef this,
                                  FunV_tFunctionalValue  iValue,
                                  xColor3fRef            iBaseColor,
                                  xColor3fRef            oColor );

/* returns the average value of the voxels in the overlay selection range.
   also returns the func idx and ras coords of that value if voxels are
   passed in. */
FunV_tErr FunV_GetAvgFunctionalValue ( tkmFunctionalVolumeRef this,
                                       FunV_tFunctionalValue* oValue,
                                       xVoxelRef              oFuncIdx,
                                       xVoxelRef              oFuncRAS,
                                       tBoolean*              obIsSelection );

/* converts an anatomical index to a functional overlay index or RAS */
FunV_tErr FunV_ConvertMRIIdxToFuncIdx ( tkmFunctionalVolumeRef this,
                                        xVoxelRef              iMRIIdx,
                                        xVoxelRef              oFuncIdx );

FunV_tErr FunV_ConvertMRIIdxToFuncRAS ( tkmFunctionalVolumeRef this,
                                        xVoxelRef              iMRIIdx,
                                        xVoxelRef              oFuncRAS );

/* time course graph access */

/* reads the graph tcl script and sets up the window. also sets up initial
   values for functional state. */
FunV_tErr FunV_InitGraphWindow ( tkmFunctionalVolumeRef this,
                                 Tcl_Interp*            pInterp );

/* for displaying voxels in graph */
FunV_tErr FunV_BeginSelectionRange      ( tkmFunctionalVolumeRef this );
FunV_tErr FunV_AddMRIIdxToSelectionRange( tkmFunctionalVolumeRef this,
    xVoxelRef              ipVoxel );
FunV_tErr FunV_EndSelectionRange        ( tkmFunctionalVolumeRef this );

/* finds average time course values for a condition over the voxels in the
   current selection range. also returns the number of good voxels */
FunV_tErr FunV_CalcTimeCourseAverages_ ( tkmFunctionalVolumeRef this,
    int                    inCondition,
    int*                   onNumValues,
    float*                 oafValues,
    float*                 oafDeviations);
/* prints out a log of the selected time course data */
FunV_tErr FunV_PrintSelectionRangeToFile ( tkmFunctionalVolumeRef this,
    char*                  isFileName );

/* grabs values for the current selected voxels and shoots them towards
   the graph to be drawn onto the screen. */
FunV_tErr FunV_DrawGraph           ( tkmFunctionalVolumeRef this );
FunV_tErr FunV_SendGraphData_      ( tkmFunctionalVolumeRef this,
                                     int                    inCondition,
                                     int                    inNumValues,
                                     float*                 iafValues );
FunV_tErr FunV_SendGraphErrorBars_ ( tkmFunctionalVolumeRef this,
                                     int                    inCondition,
                                     int                    inNumValues,
                                     float*                 iafBars );

/* manually set the location string in the graph */
FunV_tErr FunV_SetLocationString ( tkmFunctionalVolumeRef this,
                                   char*                  isLabel );

/* gets value at a point and then selects all voxels around it with a
   value >= to the starting value. */
typedef struct {
  tkmFunctionalVolumeRef mThis;
  FunV_tFunctionalValue  mStartValue;
  FunV_tFindStatsComp    mCompareType;
  int                    mnCount;
}
FunV_tFloodSelectCallbackData;

FunV_tErr FunV_FloodSelect ( tkmFunctionalVolumeRef this,
                             xVoxelRef              iSeedMRIIdx,
                             tkm_tVolumeType        iVolume,
                             int                    inDistance,
                             FunV_tFindStatsComp    iCompare );

tBoolean FunV_CompareMRIAndFuncValues ( xVoxelRef  iMRIIdx,
                                        float      iValue,
                                        void*      ipData );

Volm_tVisitCommand FunV_FloodSelectCallback ( xVoxelRef iMRIIdx,
    float     iValue,
    void*     iData );


/* tcl commands */
int FunV_TclOlSaveRegistration    ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetRegistrationToIdentity ( ClientData iClientData,
    Tcl_Interp *ipInterp,
    int argc, char *argv[] );
int FunV_TclOlRestoreRegistration ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetTimePoint        ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetCondition        ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetDisplayFlag      ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetThreshold        ( ClientData inClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclOlSetSampleType       ( ClientData inClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclTCSetNumPreStimPoints ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclTCSetTimeResolution   ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclTCSetDisplayFlag      ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );
int FunV_TclTCPrintSelectionRangeToFile ( ClientData iClientData,
    Tcl_Interp *ipInterp,
    int argc, char *argv[] );
int FunV_TclTCPrintTimeCourseData ( ClientData iClientData,
                                    Tcl_Interp *ipInterp,
                                    int argc, char *argv[] );

/* misc */
FunV_tErr FunV_SendViewStateToTcl  ( tkmFunctionalVolumeRef this );
FunV_tErr FunV_RegisterTclCommands ( tkmFunctionalVolumeRef this,
                                     Tcl_Interp*            pInterp );
FunV_tErr FunV_GetThresholdMax     ( tkmFunctionalVolumeRef this,
                                     FunV_tFunctionalValue* oValue );
FunV_tErr FunV_GetThreshold        ( tkmFunctionalVolumeRef this,
                                     FunV_tFunctionalValue* oMin,
                                     FunV_tFunctionalValue* oMid,
                                     FunV_tFunctionalValue* oSlope );

/* these check to see if we have valid callback functions and if so,
   calls them. */
FunV_tErr FunV_OverlayChanged_        ( tkmFunctionalVolumeRef this );
FunV_tErr FunV_SendTkmeditTclCommand_ ( tkmFunctionalVolumeRef this,
                                        tkm_tTclCommand        iCommand,
                                        char*                  isArguments );
FunV_tErr FunV_SendTclCommand_        ( tkmFunctionalVolumeRef this,
                                        FunV_tTclCommand       iCommand,
                                        char*                  isArguments );

/* error checking etc */
FunV_tErr FunV_Verify ( tkmFunctionalVolumeRef this );
void FunV_Signal ( char* isFuncName, int inLineNum, FunV_tErr ieCode );
char * FunV_GetErrorString ( FunV_tErr inErr );

/* dng */
void FunV_SetBlendActivation(tBoolean bFlag);

#endif
