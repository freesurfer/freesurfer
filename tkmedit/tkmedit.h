#ifndef tkmedit_h
#define tkmedit_h

//#define XDEBUG_NO_CODE

#include "mrisurf.h" /* declares WM_MIN_VAL */
#include "mriTypes.h"
#include "xVoxel.h"
#include "x3DList.h"
#include "xGLutWindow.h"
#include "mriHeadPointList.h"
#include "mriSurface.h"

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
  tkm_tErr_CouldntLoadSurfaceAnnotation,
  tkm_tErr_CouldntLoadColorTable,
  tkm_tErr_CouldntLoadHeadPointsList,
  tkm_tErr_CouldntLoadSegmentation,
  tkm_tErr_CouldntLoadOverlay,
  tkm_tErr_CouldntLoadTimeCourse,
  tkm_tErr_CouldntLoadTransform,
  tkm_tErr_CouldntLoadGCA,
  tkm_tErr_CouldntLoadVLI,
  tkm_tErr_CouldntLoadDTIVolume,
  tkm_tErr_CouldntLoadVectorField,
  tkm_tErr_ErrorAccessingFile,
  tkm_tErr_ErrorAccessingVolume,
  tkm_tErr_ErrorAccessingSegmentationVolume,
  tkm_tErr_ErrorAccessingFunctionalVolume,
  tkm_tErr_ErrorAccessingList,
  tkm_tErr_ErrorAccessingSurface,
  tkm_tErr_CouldntWriteFile,
  tkm_tErr_CouldntAllocate,
  tkm_tErr_AnatomicalVolumeNotLoaded,
  tkm_tErr_SurfaceNotLoaded,
  tkm_tErr_OverlayNotLoaded,
  tkm_tErr_GCANotLoaded,
  tkm_tErr_SegmentationNotLoaded,
  tkm_tErr_DTIVolumesDifferentSize,
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
  tkm_tTclCommand_UpdateSegLabel,
  tkm_tTclCommand_UpdateAuxSegLabel,
  tkm_tTclCommand_UpdateHeadPointLabel,
  tkm_tTclCommand_UpdateDistance,
  tkm_tTclCommand_UpdateZoomLevel,
  tkm_tTclCommand_UpdateOrientation,
  tkm_tTclCommand_UpdateDisplayFlag,
  tkm_tTclCommand_UpdateTool,
  tkm_tTclCommand_UpdateBrushTarget,
  tkm_tTclCommand_UpdateBrushShape,
  tkm_tTclCommand_UpdateBrushInfo,
  tkm_tTclCommand_UpdateFloodSelectParams,
  tkm_tTclCommand_UpdateCursorColor,
  tkm_tTclCommand_UpdateCursorShape,
  tkm_tTclCommand_UpdateSurfaceLineWidth,
  tkm_tTclCommand_UpdateSurfaceLineColor,
  tkm_tTclCommand_UpdateSegBrushInfo,
  tkm_tTclCommand_UpdateVolumeColorScale,
  tkm_tTclCommand_UpdateSegmentationVolumeAlpha,
  tkm_tTclCommand_UpdateDTIVolumeAlpha,
  tkm_tTclCommand_UpdateTimerStatus,
  tkm_tTclCommand_UpdateSubjectDirectory,
  tkm_tTclCommand_UpdateSegmentationColorTable,
  tkm_tTclCommand_UpdateVolumeDirty,
  tkm_tTclCommand_UpdateAuxVolumeDirty,
  tkm_tTclCommand_UpdateVolumeValueMinMax,
  tkm_tTclCommand_UpdateVolumeSampleType,
  tkm_tTclCommand_UpdateVolumeResampleMethod,
  
  /* display status */
  tkm_tTclCommand_ShowVolumeCoords,
  tkm_tTclCommand_ShowRASCoords,
  tkm_tTclCommand_ShowTalCoords,
  tkm_tTclCommand_ShowAuxValue,
  tkm_tTclCommand_ShowSegLabel,
  tkm_tTclCommand_ShowAuxSegLabel,
  tkm_tTclCommand_ShowHeadPointLabel,
  tkm_tTclCommand_ShowFuncCoords,
  tkm_tTclCommand_ShowFuncValue,
  tkm_tTclCommand_ShowAuxVolumeOptions,
  tkm_tTclCommand_ShowVolumeDirtyOptions,
  tkm_tTclCommand_ShowAuxVolumeDirtyOptions,
  tkm_tTclCommand_ShowMainTransformLoadedOptions,
  tkm_tTclCommand_ShowAuxTransformLoadedOptions,
  tkm_tTclCommand_ShowFuncOverlayOptions,
  tkm_tTclCommand_ShowFuncTimeCourseOptions,
  tkm_tTclCommand_ShowSurfaceLoadingOptions,
  tkm_tTclCommand_ShowSurfaceViewingOptions,
  tkm_tTclCommand_ShowOriginalSurfaceViewingOptions,
  tkm_tTclCommand_ShowPialSurfaceViewingOptions,
  tkm_tTclCommand_ShowHeadPointLabelEditingOptions,
  tkm_tTclCommand_ShowVLIOptions,
  tkm_tTclCommand_ShowGCAOptions,
  tkm_tTclCommand_ShowDTIOptions,
  tkm_tTclCommand_ShowOverlayRegistrationOptions,
  tkm_tTclCommand_ShowSegmentationOptions,
  tkm_tTclCommand_ShowAuxSegmentationOptions,
  tkm_tTclCommand_ClearSegColorTable,
  tkm_tTclCommand_AddSegColorTableEntry,
  
  /* histogram */
  tkm_tTclCommand_DrawHistogram,
  
  /* interface configuration */
  tkm_tTclCommand_MoveToolWindow,
  tkm_tTclCommand_RaiseToolWindow,
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
  tkm_knNumVolumeTypes
} tkm_tVolumeType;

typedef enum {
  tkm_tSegType_Main = 0,
  tkm_tSegType_Aux,
  tkm_knNumSegTypes
} tkm_tSegType;

typedef enum {
  tkm_tVolumeTarget_MainAna = 0,
  tkm_tVolumeTarget_AuxAna,
  tkm_tVolumeTarget_MainSeg,
  tkm_tVolumeTarget_AuxSeg,
  tkm_knNumVolumeTargets
} tkm_tVolumeTarget;

typedef enum {
  tkm_tSurfaceType_Main = 0,
  tkm_tSurfaceType_Aux,
  tkm_knNumSurfaceTypes
} tkm_tSurfaceType;

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
void tkm_GetAnaDimension  ( tkm_tVolumeType iVolume,
			    int*            onDimensionX, 
			    int*            onDimensionY, 
			    int*            onDimensionZ );
tBoolean tkm_IsValidAnaIdx ( tkm_tVolumeType iVolume,
			     xVoxelRef       iAnaIdx );

/* segmentation value */
void tkm_GetSegmentationColorAtVoxel ( tkm_tSegType iVolume,
				       xVoxelRef    iAnaIdx,
				       xColor3fRef  iBaseColor,
				       xColor3fRef  oColor );
void tkm_GetSegLabel                 ( tkm_tSegType iVolume,
				       xVoxelRef    iAnaIdx, 
				       int*         onIndex,
				       char*        osLabel );

/* selects all the voxels in the label with the given index */
void tkm_SelectCurrentSegLabel     ( tkm_tSegType iVolume,
				     int          inIndex );
/* graphs the avg of all the voxels in the label with the given index */
void tkm_GraphCurrentSegLabelAvg   ( tkm_tSegType iVolume,
				     int          inIndex );

/* get the volume of an segmentation label */
void tkm_CalcSegLabelVolume ( tkm_tSegType iVolume,
			      xVoxelRef    iAnaIdx,
			      int*         onVolume );
 
/* editing the segmentation */
void tkm_EditSegmentation      ( tkm_tSegType      iVolume,
				 xVoxelRef         iMRIIdx,
				 int               inIndex );
void tkm_EditSegmentationArray ( tkm_tSegType      iVolume,
				 xVoxelRef         iaMRIIdx,
				 int               inCount,
				 int               inIndex );
void tkm_FloodFillSegmentation ( tkm_tSegType      iVolume,
				 xVoxelRef         iMRIIdx,
				 int               inIndex,
				 tBoolean          ib3D,
				 tkm_tVolumeTarget iSrc,
				 int               inFuzzy,
				 int               inDistance );

/* dti color */
void tkm_GetDTIColorAtVoxel ( xVoxelRef        iAnaIdx,
			      mri_tOrientation iPlane,
			      xColor3fRef      iBaseColor,
			      xColor3fRef      oColor );
            

/* dealing with control points */
void tkm_MakeControlPoint             ( xVoxelRef        iAnaIdx );
void tkm_RemoveControlPointWithinDist ( xVoxelRef        iAnaIdx,
					mri_tOrientation iPlane,
					int              inDistance );
void tkm_WriteControlFile             ();

/* editing */
void tkm_EditAnatomicalVolumeInRange( tkm_tVolumeType  iVolume, 
				      xVoxelRef        iMRIIdx, 
				      tVolumeValue     inLow, 
				      tVolumeValue     inHigh, 
				      tVolumeValue     inNewValue );
void tkm_EditAnatomicalVolumeInRangeArray( tkm_tVolumeType  iVolume, 
					   xVoxelRef        iaMRIIdx, 
					   int              inCount,
					   tVolumeValue     inLow, 
					   tVolumeValue     inHigh, 
					   tVolumeValue     inNewValue );

/* Sets a region in the anatomical volume to a new value. */
void tkm_SetAnatomicalVolumeRegion ( tkm_tVolumeType iVolume,
				     int             iAnaX0,
				     int             iAnaX1,
				     int             iAnaY0,
				     int             iAnaY1,
				     int             iAnaZ0,
				     int             iAnaZ1,
				     float           iNewValue );

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
void tkm_SelectVoxel         ( xVoxelRef iMRIIdx );
void tkm_SelectVoxelArray    ( xVoxelRef iaMRIIdx, int inCount );
void tkm_DeselectVoxel       ( xVoxelRef iMRIIdx );
void tkm_DeselectVoxelArray  ( xVoxelRef iaMRIIdx, int inCount );
void tkm_ClearSelection      ();
tBoolean tkm_IsSelectionPresent ();
void tkm_FloodSelect         ( xVoxelRef         iSeedAnaIdx,
			       tBoolean          ib3D,
			       tkm_tVolumeTarget iSrc,
			       int               inFuzzy,
			       int               inDistance,
			       tBoolean          ibSelect );

/* event processing */
void tkm_HandleIdle ();

/* writing points out to files. */
void tkm_WriteVoxelToControlFile ( xVoxelRef iAnaIdx );
void tkm_WriteVoxelToEditFile    ( xVoxelRef iAnaIdx );
void tkm_ReadCursorFromEditFile  ();

/* writing surface distances. */
void tkm_SetSurfaceDistance    ( xVoxelRef iAnaIdx,
				 float     ifDistance );

/* show nearest verts */
void tkm_ShowNearestSurfaceVertex ( Surf_tVertexSet iVertexSet );

/* cleaning up */
void tkm_Quit ();

/* send a tcl command */
void tkm_SendTclCommand ( tkm_tTclCommand iCommand,
        char*           isArguments );

char* tkm_GetErrorString( tkm_tErr ieCode );

#endif

