#ifndef tkmedit_h
#define tkmedit_h

#include "tkmVoxel.h"
#include "tkmVoxelSpace.h"
#include "tkmVoxelList.h"
#include "mrisurf.h"
#include "xGLutWindow.h"

typedef unsigned char tVolumeValue;
typedef tVolumeValue* tVolumeRef;

#define knMinVolumeValue 0
#define knMaxVolumeValue 255

#define WM_EDITED_OFF 1
#define WM_MIN_VAL    2 /* 1 is used for voxels that are edited to off */

#define tkm_knEditToWhiteLow      knMinVolumeValue
#define tkm_knEditToWhiteHigh     WM_MIN_VAL
#define tkm_knEditToWhiteNewValue 255

#define tkm_knEditToBlackLow      knMinVolumeValue
#define tkm_knEditToBlackHigh     knMaxVolumeValue
#define tkm_knEditToBlackNewValue WM_EDITED_OFF

/* this is synced to defs in tkmedit.c and in tkm_interface.tcl */
typedef enum {  
  tkm_tSurfaceType_None = -1,
  tkm_tSurfaceType_Current = 0,
  tkm_tSurfaceType_Original,
  tkm_tSurfaceType_Canonical,
  tkm_knNumSurfaceTypes
} tkm_tSurfaceType;

/* this is synced to defs in tkmedit.c and in tkm_interface.tcl */
typedef enum { 
  tkm_tOrientation_None = -1,
  tkm_tOrientation_Coronal = 0,
  tkm_tOrientation_Horizontal,
  tkm_tOrientation_Sagittal,
  tkm_knNumOrientations
} tkm_tOrientation;

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

  /* display status */
  tkm_tTclCommand_ShowVolumeCoords,
  tkm_tTclCommand_ShowRASCoords,
  tkm_tTclCommand_ShowTalCoords,
  tkm_tTclCommand_ShowAuxValue,
  tkm_tTclCommand_ShowFuncCoords,
  tkm_tTclCommand_ShowFuncValue,
  tkm_tTclCommand_ShowFuncOverlayOptions,
  tkm_tTclCommand_ShowFuncTimeCourseOptions,

  /* interface configuration */
  tkm_tTclCommand_MoveToolWindow,
  tkm_tTclCommand_CsurfInterface,
  tkm_tTclCommand_ErrorDlog,
  tkm_knNumTclCommands
} tkm_tTclCommand;

/* convesion functions */
void tkm_ConvertVolumeToRAS ( VoxelRef inVolumeVox, VoxelRef outRASVox );
void tkm_ConvertVolumeToTal ( VoxelRef inVolumeVox, VoxelRef outTalVox );

/* interfaces for accessing current volume state */
unsigned char tkm_GetVolumeValue ( tVolumeRef, VoxelRef );
void tkm_GetAnatomicalVolumeColor( tVolumeValue inValue,
           float* outRed,
           float* outGreen, 
           float* outBlue );

/* parcellation value */
void tkm_GetParcellationColor( VoxelRef, 
             float* outRed,
             float* outGreen, 
             float* outBlue );

/* dealing with control points */
void tkm_AddNearestCtrlPtToSelection ( VoxelRef inVolumeVox, 
               tkm_tOrientation inPlane );
void tkm_RemoveNearestCtrlPtFromSelection ( VoxelRef inVolumeVox, 
              tkm_tOrientation inPlane );
void tkm_NewCtrlPt ();
void tkm_DeselectAllCtrlPts ();
void tkm_DeleteSelectedCtrlPts ();
void tkm_WriteControlFile ();

/* editing */
void tkm_EditVoxelInRange( VoxelRef     inVolumeVox, 
         tVolumeValue inLow, 
         tVolumeValue inHigh, 
         tVolumeValue inNewValue );

/* undo list */
void tkm_ClearUndoList ();
void tkm_RestoreUndoList ();

/* selecting */
void tkm_SelectVoxel ( VoxelRef inVolumeVox );
void tkm_DeselectVoxel ( VoxelRef inVolumeVox );
void tkm_ClearSelection ();

/* event processing */
void tkm_HandleIdle ();

/* writing points out to files. */
void tkm_WriteVoxelToControlFile ( VoxelRef inVolumeVox );
void tkm_WriteVoxelToEditFile ( VoxelRef inVolumeVox );
void tkm_ReadCursorFromEditFile ();

/* cleaning up */
void tkm_Quit ();

/* various global variable access */
char* tkm_GetSubjectName();
char* tkm_GetVolumeName();
char* tkm_GetAuxVolumeName();

/* send a tcl command */
void tkm_SendTclCommand ( tkm_tTclCommand inCommand,
        char* inArguments );

#endif
