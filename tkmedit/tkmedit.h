#ifndef tkmedit_h
#define tkmedit_h

#include "mriTypes.h"
#include "xVoxel.h"
#include "x3DList.h"
#include "mrisurf.h"
#include "xGLutWindow.h"

typedef unsigned char tVolumeValue;
typedef tVolumeValue* tVolumeRef;

#define knMinVolumeValue 0
#define knMaxVolumeValue 255

#define WM_EDITED_OFF 1

#ifndef WM_MIN_VAL
#define WM_MIN_VAL    5 /* 1 is used for voxels that are edited to off */
#endif

#define tkm_knEditToWhiteLow      knMinVolumeValue
#define tkm_knEditToWhiteHigh     WM_MIN_VAL
#define tkm_knEditToWhiteNewValue 255

#define tkm_knEditToBlackLow      knMinVolumeValue
#define tkm_knEditToBlackHigh     knMaxVolumeValue
#define tkm_knEditToBlackNewValue WM_EDITED_OFF

/* this is synced to defs in tkmedit.c and in tkm_interface.tcl */
/*
typedef enum {  
  Surf_tVertexSet_None = -1,
  Surf_tVertexSet_Main = 0,
  Surf_tVertexSet_Original,
  Surf_tVertexSet_Pial,
  Surf_knNumVertexSets
} Surf_tVertexSet;
*/

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
  tkm_tTclCommand_UpdateBrush,
  tkm_tTclCommand_UpdateBrushThreshold,
  tkm_tTclCommand_UpdateVolumeColorScale,
  tkm_tTclCommand_UpdateParcellationLabel,

  /* display status */
  tkm_tTclCommand_ShowVolumeCoords,
  tkm_tTclCommand_ShowRASCoords,
  tkm_tTclCommand_ShowTalCoords,
  tkm_tTclCommand_ShowAuxValue,
  tkm_tTclCommand_ShowParcellationLabel,
  tkm_tTclCommand_ShowFuncCoords,
  tkm_tTclCommand_ShowFuncValue,
  tkm_tTclCommand_ShowFuncOverlayOptions,
  tkm_tTclCommand_ShowFuncTimeCourseOptions,
  tkm_tTclCommand_ShowSurfaceLoadingOptions,
  tkm_tTclCommand_ShowOriginalSurfaceViewingOptions,
  tkm_tTclCommand_ShowCanonicalSurfaceViewingOptions,

  /* interface configuration */
  tkm_tTclCommand_MoveToolWindow,
  tkm_tTclCommand_CsurfInterface,
  tkm_tTclCommand_ErrorDlog,
  tkm_knNumTclCommands
} tkm_tTclCommand;

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
void tkm_GetAnatomicalVolumeColor( tVolumeValue inValue, 
           xColor3fRef  oColor );

/* getting the maximum intensity projection */
tVolumeValue tkm_GetMaxIntProjValue( tVolumeRef       iVolume, 
             mri_tOrientation iOrientation, 
             xVoxelRef        ipVoxel );

/* parcellation value */
void tkm_GetParcellationColor ( xVoxelRef   iWhere, 
        xColor3fRef oColor );
void tkm_GetParcellationLabel ( xVoxelRef   iWhere, 
        char*       osLabel );

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
