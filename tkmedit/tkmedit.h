#ifndef tkmedit_h
#define tkmedit_h

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

/* commands for the tcl side of things */
typedef enum {

  /* updating vars */
  tkm_tTclCommand_UpdateLinkedCursorFlag = 0,
  tkm_tTclCommand_UpdateVolumeCursor,
  tkm_tTclCommand_UpdateRASCursor,
  tkm_tTclCommand_UpdateTalCursor,
  tkm_tTclCommand_UpdateVolumeName,
  tkm_tTclCommand_UpdateVolumeValue,
  tkm_tTclCommand_UpdateAuxVolumeName,
  tkm_tTclCommand_UpdateAuxVolumeValue,
  tkm_tTclCommand_UpdateFunctionalCoords,
  tkm_tTclCommand_UpdateFunctionalValue,
  tkm_tTclCommand_UpdateZoomLevel,
  tkm_tTclCommand_UpdateOrientation,
  tkm_tTclCommand_UpdateDisplayFlag,
  tkm_tTclCommand_UpdateTool,
  tkm_tTclCommand_UpdateBrushShape,
  tkm_tTclCommand_UpdateBrushInfo,
  tkm_tTclCommand_UpdateVolumeColorScale,
  tkm_tTclCommand_UpdateROILabel,
  tkm_tTclCommand_UpdateHeadPointLabel,

  /* display status */
  tkm_tTclCommand_ShowVolumeCoords,
  tkm_tTclCommand_ShowRASCoords,
  tkm_tTclCommand_ShowTalCoords,
  tkm_tTclCommand_ShowAuxValue,
  tkm_tTclCommand_ShowROILabel,
  tkm_tTclCommand_ShowHeadPointLabel,
  tkm_tTclCommand_ShowFuncCoords,
  tkm_tTclCommand_ShowFuncValue,
  tkm_tTclCommand_ShowFuncOverlayOptions,
  tkm_tTclCommand_ShowFuncTimeCourseOptions,
  tkm_tTclCommand_ShowSurfaceLoadingOptions,
  tkm_tTclCommand_ShowSurfaceViewingOptions,
  tkm_tTclCommand_ShowOriginalSurfaceViewingOptions,
  tkm_tTclCommand_ShowCanonicalSurfaceViewingOptions,
  tkm_tTclCommand_ShowHeadPointLabelEditingOptions,
  tkm_tTclCommand_ShowOverlayRegistrationOptions,

  /* interface configuration */
  tkm_tTclCommand_MoveToolWindow,
  tkm_tTclCommand_CsurfInterface,

  /* misc */
  tkm_tTclCommand_ErrorDlog,  
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

  tkm_tAxis_X = 0,
  tkm_tAxis_Y,
  tkm_tAxis_Z,
  tkm_knNumAxes
} tkm_tAxis;

/* convesion functions */
void tkm_ConvertVolumeToRAS ( xVoxelRef iAnaIdx,
            xVoxelRef oRASVox );
void tkm_ConvertRASToVolume ( xVoxelRef iRASVox,
            xVoxelRef oAnaIdx );
void tkm_ConvertVolumeToTal ( xVoxelRef iAnaIdx,
            xVoxelRef oTalVox );
void tkm_ConvertTalToVolume ( xVoxelRef iTalVox,
            xVoxelRef oAnaIdx );

tBoolean tkm_IsValidVolumeIdx ( xVoxelRef iAnaIdx );
tBoolean tkm_IsValidRAS       ( xVoxelRef iRAS );

/* interfaces for accessing current volume state */
tVolumeValue tkm_GetVolumeValue ( tVolumeRef iVolume,
          xVoxelRef  iAnaIdx );
void tkm_GetAnatomicalVolumeColor( tVolumeRef   iVolume,
           tVolumeValue inValue, 
           xColor3fRef  oColor );

/* getting the maximum intensity projection */
tVolumeValue tkm_GetMaxIntProjValue( tVolumeRef       iVolume, 
             mri_tOrientation iOrientation, 
             xVoxelRef        ipVoxel );

/* roi value */
void tkm_GetROIColorAtVoxel ( xVoxelRef   iWhere, 
            xColor3fRef oColor );
void tkm_GetROILabel ( xVoxelRef   iWhere, 
           int*        onIndex,
           char*       osLabel );

/* selects all the voxels in the label with the given index */
void tkm_SelectCurrentROI ( int inIndex );

/* dealing with control points */
void tkm_AddNearestCtrlPtToSelection      ( xVoxelRef        iAnaIdx, 
              mri_tOrientation iPlane );
void tkm_RemoveNearestCtrlPtFromSelection ( xVoxelRef        iAnaIdx, 
              mri_tOrientation iPlane );
void tkm_NewCtrlPt                        ();
void tkm_DeselectAllCtrlPts               ();
void tkm_DeleteSelectedCtrlPts            ();
void tkm_WriteControlFile                 ();

/* editing */
void tkm_EditVoxelInRange( xVoxelRef    iAnaIdx, 
         tVolumeValue inLow, 
         tVolumeValue inHigh, 
         tVolumeValue inNewValue );

/* undo list */
void tkm_ClearUndoList   ();
void tkm_RestoreUndoList ();

/* head points */
void tkm_GetHeadPoint ( xVoxelRef           iAnaIdx,
      mri_tOrientation    iOrientation,
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

/* various global variable access */
char* tkm_GetSubjectName  ();
char* tkm_GetVolumeName   ();
char* tkm_GetAuxVolumeName();

/* send a tcl command */
void tkm_SendTclCommand ( tkm_tTclCommand iCommand,
        char*           isArguments );

#endif
