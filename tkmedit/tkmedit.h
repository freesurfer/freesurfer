#ifndef tkmedit_h
#define tkmedit_h

//#define XDEBUG_NO_CODE

#include "mrisurf.h" /* declares WM_MIN_VAL */
#include "mriTypes.h"
#include "xVoxel.h"
#include "x3DList.h"
#include "xGLutWindow.h"
#include "mriHeadPointList.h"

typedef unsigned char tVolumeValue;
typedef tVolumeValue* tVolumeRef;

#define knMinVolumeValue 0
#define knMaxVolumeValue 255

#define WM_EDITED_OFF 1

#define tkm_knEditToWhiteLow      knMinVolumeValue
#define tkm_knEditToWhiteHigh     WM_MIN_VAL
#define tkm_knEditToWhiteNewValue 255

#define tkm_knEditToBlackLow      WM_MIN_VAL
#define tkm_knEditToBlackHigh     knMaxVolumeValue
#define tkm_knEditToBlackNewValue WM_EDITED_OFF

#define tkm_knErrStringLen   256
#define tkm_knPathLen        1024
#define tkm_knNameLen        256
#define tkm_knTclCmdLen      1024

typedef enum {
  tkm_tErr_InvalidErrorCode = -1,
  tkm_tErr_NoErr,
  tkm_tErr_InvalidParameter,
  tkm_tErr_MRIDIRNotDefined,
  tkm_tErr_SUBJECTSDIRNotDefined,
  tkm_tErr_CouldntGetPWD,
  tkm_tErr_CouldntInitTcl,
  tkm_tErr_CouldntInitTk,
  tkm_tErr_CouldntInitTix,
  tkm_tErr_CouldntInitBLT,
  tkm_tErr_CouldntReadVolume,
  tkm_tErr_CouldntLoadSurface,
  tkm_tErr_CouldntLoadLabel,  
  tkm_tErr_CouldntLoadSurfaceVertexSet,
  tkm_tErr_CouldntLoadColorTable,
  tkm_tErr_CouldntLoadHeadPointsList,
  tkm_tErr_CouldntLoadParcellation,
  tkm_tErr_CouldntLoadOverlay,
  tkm_tErr_CouldntLoadTimeCourse,
  tkm_tErr_CouldntLoadTransform,
  tkm_tErr_CouldntLoadGCA,
  tkm_tErr_ErrorAccessingFile,
  tkm_tErr_ErrorAccessingVolume,
  tkm_tErr_ErrorAccessingROIGroup,
  tkm_tErr_ErrorAccessingList,
  tkm_tErr_CouldntWriteFile,
  tkm_tErr_CouldntAllocate,
  tkm_tErr_SurfaceNotLoaded,
  tkm_tErr_CouldntCacheScriptName,
  tkm_tErr_InvalidScriptName,
  tkm_tErr_GetTimeOfDayFailed,
  tkm_tErr_Unrecoverable,
  tkm_knNumErrorCodes
} tkm_tErr;

/* commands for the tcl side of things */
typedef enum {

  /* updating vars */
  tkm_tTclCommand_UpdateLinkedCursorFlag = 0,
  tkm_tTclCommand_UpdateVolumeCursor,
  tkm_tTclCommand_UpdateVolumeSlice,
  tkm_tTclCommand_UpdateRASCursor,
  tkm_tTclCommand_UpdateTalCursor,
  tkm_tTclCommand_UpdateScannerCursor,
  tkm_tTclCommand_UpdateMNICursor,
  tkm_tTclCommand_UpdateVolumeName,
  tkm_tTclCommand_UpdateVolumeValue,
  tkm_tTclCommand_UpdateAuxVolumeName,
  tkm_tTclCommand_UpdateAuxVolumeValue,
  tkm_tTclCommand_UpdateFunctionalCoords,
  tkm_tTclCommand_UpdateFunctionalRASCoords,
  tkm_tTclCommand_UpdateFunctionalValue,
  tkm_tTclCommand_UpdateZoomLevel,
  tkm_tTclCommand_UpdateOrientation,
  tkm_tTclCommand_UpdateDisplayFlag,
  tkm_tTclCommand_UpdateTool,
  tkm_tTclCommand_UpdateBrushShape,
  tkm_tTclCommand_UpdateBrushInfo,
  tkm_tTclCommand_UpdateCursorColor,
  tkm_tTclCommand_UpdateCursorShape,
  tkm_tTclCommand_UpdateSurfaceLineWidth,
  tkm_tTclCommand_UpdateSurfaceLineColor,
  tkm_tTclCommand_UpdateParcBrushInfo,
  tkm_tTclCommand_UpdateVolumeColorScale,
  tkm_tTclCommand_UpdateROILabel,
  tkm_tTclCommand_UpdateROIGroupAlpha,
  tkm_tTclCommand_UpdateHeadPointLabel,
  tkm_tTclCommand_UpdateTimerStatus,
  tkm_tTclCommand_UpdateHomeDirectory,

  /* display status */
  tkm_tTclCommand_ShowVolumeCoords,
  tkm_tTclCommand_ShowRASCoords,
  tkm_tTclCommand_ShowTalCoords,
  tkm_tTclCommand_ShowAuxValue,
  tkm_tTclCommand_ShowROILabel,
  tkm_tTclCommand_ShowHeadPointLabel,
  tkm_tTclCommand_ShowFuncCoords,
  tkm_tTclCommand_ShowFuncValue,
  tkm_tTclCommand_ShowAuxVolumeOptions,
  tkm_tTclCommand_ShowFuncOverlayOptions,
  tkm_tTclCommand_ShowFuncTimeCourseOptions,
  tkm_tTclCommand_ShowSurfaceLoadingOptions,
  tkm_tTclCommand_ShowSurfaceViewingOptions,
  tkm_tTclCommand_ShowOriginalSurfaceViewingOptions,
  tkm_tTclCommand_ShowCanonicalSurfaceViewingOptions,
  tkm_tTclCommand_ShowHeadPointLabelEditingOptions,
  tkm_tTclCommand_ShowROIGroupOptions,
  tkm_tTclCommand_ShowOverlayRegistrationOptions,
  tkm_tTclCommand_ShowParcellationOptions,
  tkm_tTclCommand_ClearParcColorTable,
  tkm_tTclCommand_AddParcColorTableEntry,

  /* interface configuration */
  tkm_tTclCommand_MoveToolWindow,
  tkm_tTclCommand_CsurfInterface,
  tkm_tTclCommand_FinishBuildingInterface,

  /* misc */
  tkm_tTclCommand_ErrorDlog,  
  tkm_tTclCommand_FormattedErrorDlog,
  tkm_tTclCommand_AlertDlog,
  tkm_tTclCommand_MakeProgressDlog,
  tkm_tTclCommand_UpdateProgressDlog,
  tkm_tTclCommand_DestroyProgressDlog,
  tkm_knNumTclCommands
} tkm_tTclCommand;

typedef enum {

  tkm_tVolumeType_Main = 0,
  tkm_tVolumeType_Aux,
  tkm_tVolumeType_Parc,
  tkm_knNumVolumeTypes
} tkm_tVolumeType;

typedef enum {

  tkm_tAxis_X = 0,
  tkm_tAxis_Y,
  tkm_tAxis_Z,
  tkm_knNumAxes
} tkm_tAxis;

// ==================================================================== OUTPUT

#define InitOutput
#define DeleteOutput
#define OutputPrint            fprintf ( stdout,
#define EndOutputPrint         ); fflush( stdout );

// ===========================================================================

/* progress bar functions */
void tkm_MakeProgressBar   ( char* isName, char* isDesc );
void tkm_UpdateProgressBar ( char* isName, float ifPercent );
void tkm_FinishProgressBar ( char* isName );

/* output functions */
void tkm_DisplayMessage ( char* isMessage );
void tkm_DisplayError   ( char* isAction, char* isError, char* isDesc );
void tkm_DisplayAlert   ( char* isAction, char* isMsg, char* isDesc );

/* volume value */
void tkm_GetValueAtAnaIdx ( tkm_tVolumeType iVolume,
          xVoxelRef       iAnaIdx,
          tVolumeValue*   oValue );
void tkm_GetAnaDimension  ( tkm_tVolumeType iVolume,
          int*            onDimension );
tBoolean tkm_IsValidAnaIdx ( tkm_tVolumeType iVolume,
           xVoxelRef       iAnaIdx );

/* roi value */
void tkm_GetROIColorAtVoxel ( xVoxelRef   iAnaIdx,
            xColor3fRef iBaseColor,
            xColor3fRef oColor );
void tkm_GetROILabel        ( xVoxelRef   iAnaIdx, 
            int*        onIndex,
            char*       osLabel );

/* selects all the voxels in the label with the given index */
void tkm_SelectCurrentROI   ( int inIndex );

/* get the volume of an roi */
void tkm_CalcROIVolume ( xVoxelRef iAnaIdx,
       int*      onVolume );

/* editing the parcellation */
void tkm_EditParcellation      ( xVoxelRef       iAnaIdx,
         int             inIndex );
void tkm_FloodFillParcellation ( xVoxelRef       iAnaIdx,
         int             inIndex,
         tBoolean        ib3D,
         tkm_tVolumeType iSrc,
         int             inFuzzy,
         int             inDistance );

/* dealing with control points */
void tkm_MakeControlPoint             ( xVoxelRef        iAnaIdx );
void tkm_RemoveControlPointWithinDist ( xVoxelRef        iAnaIdx,
          mri_tOrientation iPlane,
          int              inDistance );
void tkm_WriteControlFile             ();

/* editing */
void tkm_EditVoxelInRange( xVoxelRef    iAnaIdx, 
         tVolumeValue inLow, 
         tVolumeValue inHigh, 
         tVolumeValue inNewValue );

/* undo list */
void tkm_ClearUndoList   ();
void tkm_RestoreUndoList ();

/* undo volume */
void     tkm_ClearUndoVolume               ();
void     tkm_RestoreUndoVolumeAroundAnaIdx ( xVoxelRef iAnaIdx );
tBoolean tkm_IsAnaIdxInUndoVolume          ( xVoxelRef iAnaIdx );

/* head points */
void tkm_GetHeadPoint ( xVoxelRef           iAnaIdx,
      mri_tOrientation    iOrientation,
      tBoolean            ibFlat,
      HPtL_tHeadPointRef* opPoint );

/* selecting */
void tkm_SelectVoxel    ( xVoxelRef iAnaIdx );
void tkm_DeselectVoxel  ( xVoxelRef iAnaIdx );
void tkm_ClearSelection ();

/* event processing */
void tkm_HandleIdle ();

/* writing points out to files. */
void tkm_WriteVoxelToControlFile ( xVoxelRef iAnaIdx );
void tkm_WriteVoxelToEditFile    ( xVoxelRef iAnaIdx );
void tkm_ReadCursorFromEditFile  ();

/* cleaning up */
void tkm_Quit ();

/* send a tcl command */
void tkm_SendTclCommand ( tkm_tTclCommand iCommand,
        char*           isArguments );

char* tkm_GetErrorString( tkm_tErr ieCode );

#endif

