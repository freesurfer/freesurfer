#ifndef tkmDisplayArea_h
#define tkmDisplayArea_h

#include "xTypes.h"
#include <GL/glut.h>
#include "tkmVoxel.h"
#include "tkmedit.h"
#include "xGLutWindow.h"
#include "xList.h"
#include "tkmFunctionalVolume.h"
#include "mrishash.h"

typedef enum {

  DspA_tErr_NoErr = 0,
  DspA_tErr_AllocationFailed,
  DspA_tErr_InvalidPtr,
  DspA_tErr_InvalidSignature,
  DspA_tErr_InvalidParameter,
  DspA_tErr_InvalidCursor,
  DspA_tErr_InvalidDisplayFlag,
  DspA_tErr_InvalidOrientation,
  DspA_tErr_InvalidVolumeVoxel,
  DspA_tErr_InvalidBufferPoint,
  DspA_tErr_InvalidScreenPoint,
  DspA_tErr_ErrorCreatingSurfaceHashTable,
  DspA_tErr_OutOfMemory,
  DspA_tErr_ErrorAccessingControlPoints,
  DspA_tErr_ErrorAccessingSelection,
  DspA_tErr_ErrorAccessingWindow,
  DspA_tErr_ErrorAccessingFunctionalVolume,
  DspA_tErr_ErrorAccessingSurfaceList,
  DspA_tErr_InvalidErrorCode,
  DspA_knNumErrorCodes

} DspA_tErr;

/* synced to tkm_interface.tcl */
typedef enum {

  DspA_tDisplayFlag_None = 0, 
  DspA_tDisplayFlag_AuxVolume,
  DspA_tDisplayFlag_Cursor,
  DspA_tDisplayFlag_MainSurface,
  DspA_tDisplayFlag_OriginalSurface,
  DspA_tDisplayFlag_CanonicalSurface,
  DspA_tDisplayFlag_InterpolateSurfaceVertices,
  DspA_tDisplayFlag_DisplaySurfaceVertices,
  DspA_tDisplayFlag_ControlPoints,
  DspA_tDisplayFlag_Selection,
  DspA_tDisplayFlag_FunctionalOverlay,
  DspA_tDisplayFlag_ParcellationOverlay,
  DspA_tDisplayFlag_FocusFrame,
  DspA_tDisplayFlag_Axes,
  DspA_tDisplayFlag_MaxIntProj,
  DspA_knNumDisplayFlags
  
} DspA_tDisplayFlag;

/* synced to tkm_interface.tcl */
typedef enum {

  DspA_tTool_None         = -1,
  DspA_tTool_SelectVoxels =  0,
  DspA_tTool_EditVoxels,
  DspA_tTool_SelectCtrlPts,
  DspA_tTool_CustomEditVoxels,
  DspA_knNumTools

} DspA_tTool;

/* synced to tkm_interface.tcl */
typedef enum {

  DspA_tBrushShape_None   = -1,
  DspA_tBrushShape_Circle =  0,
  DspA_tBrushShape_Square,
  DspA_knNumBrushShapes

} DspA_tBrushShape;

typedef enum {

  DspA_tVolumeType_Main = 0,
  DspA_tVolumeType_Aux,
  DspA_knNumVolumeTypes
} DspA_tVolumeType;

/* for stashing tons of surface draw lists. */
#define DspA_knMaxPointsPerPointListNode 1000
typedef struct {

  int      mnNumPoints;
  float    mafPoints[DspA_knMaxPointsPerPointListNode][2];

} tkmPointListNode, *tkmPointListNodeRef;

#define DspA_kSignature 0x194ffb2

struct tkmDisplayArea {

  tSignature        mSignature;
  
  /* superpane info */
  struct tkmMeditWindow* mpWindow;
  
  /* our size and location */
  xPoint2n          mLocationInSuper;
  int               mnWidth;
  int               mnHeight;
  
  /* frame buffer */
  int               mnVolumeSize;
  GLubyte*          mpFrameBuffer;
  float             mfFrameBufferScaleX;
  float             mfFrameBufferScaleY;

  /* view state */
  VoxelRef          mpCursor;
  tkm_tOrientation  mOrientation;
  int               mnZoomLevel;
  VoxelRef          mpZoomCenter;
  int               mnHilitedVertexIndex;
  tkm_tSurfaceType  mHilitedSurface;
  tBoolean          mabDisplayFlags [DspA_knNumDisplayFlags];
  tBoolean          mbSliceChanged;

  /* surface lists */
  xListRef*         maSurfaceLists;

  /* display data */
  tVolumeRef             mpVolume;
  tVolumeRef             mpAuxVolume;
  tVolumeRef             mpParcellationVolume;
  MRI_SURFACE*           mpSurface;
  tkmFunctionalVolumeRef mpFunctionalVolume;
  VoxelSpaceRef          mpControlPoints;
  VoxelListRef           mpSelectedControlPoints;
  VoxelSpaceRef          mpSelection;

};
typedef struct tkmDisplayArea tkmDisplayArea;
typedef tkmDisplayArea *tkmDisplayAreaRef;

#define DspA_knNumVolumeValues       256
#define DspA_knDefaultVolumeMidLevel 160;

#define DspA_knNumBytesPerPixel     4
#define DspA_knMaxPixelValue        (GLubyte)255
#define DspA_knRedPixelCompIndex    0
#define DspA_knGreenPixelCompIndex  1
#define DspA_knBluePixelCompIndex   2
#define DspA_knAlphaPixelCompIndex  3

#define DspA_knMaxZoomLevel 16
#define DspA_knMinZoomLevel  1
  
#define DspA_knMaxBrushRadius 100
#define DspA_knMinBrushRadius  1
  
#define DspA_knCursorCrosshairSize       4
#define DspA_knControlPointCrosshairSize 4
#define DspA_knSurfaceVertexSize         2

#define DspA_knSelectionIntensityIncrement 100
#define DspA_knMinSelectionIntensityDiff    10

#define DspA_kfHashTableResolution 2.0

DspA_tErr DspA_New    ( tkmDisplayAreaRef* oppWindow,
      struct tkmMeditWindow*  ipWindow );
DspA_tErr DspA_Delete ( tkmDisplayAreaRef* ioppWindow );

/* size and location in super */
DspA_tErr DspA_SetPosition ( tkmDisplayAreaRef this,
           xPoint2n          iLocation,
           int               inWidth,
           int               inHeight );

/* update window title */
DspA_tErr DspA_UpdateWindowTitle ( tkmDisplayAreaRef this );

/* all the things it will draw */
DspA_tErr DspA_SetVolume                     ( tkmDisplayAreaRef this,
                 tVolumeRef        ipVolume,
                 int               inSize );
DspA_tErr DspA_SetAuxVolume                  ( tkmDisplayAreaRef this,
                 tVolumeRef        ipVolume,
                 int               inSize );
DspA_tErr DspA_SetParcellationVolume         ( tkmDisplayAreaRef this,
                 tVolumeRef        ipVolume,
                 int               inSize );
DspA_tErr DspA_SetSurface                    ( tkmDisplayAreaRef this, 
                 MRI_SURFACE*      ipSurface );
DspA_tErr DspA_SetOverlayVolume              ( tkmDisplayAreaRef this,
                 tkmFunctionalVolumeRef ipVol );
DspA_tErr DspA_SetControlPointsSpace         ( tkmDisplayAreaRef this,
                 VoxelSpaceRef     ipVoxels );
DspA_tErr DspA_SetControlPointsSelectionList ( tkmDisplayAreaRef this,
                 VoxelListRef      ipVoxels );
DspA_tErr DspA_SetSelectionSpace             ( tkmDisplayAreaRef this, 
                 VoxelSpaceRef     ipVoxels );

/* viewing state changes */
DspA_tErr DspA_SetCursor             ( tkmDisplayAreaRef this, 
               VoxelRef          ipCursor );
DspA_tErr DspA_SetOrientation        ( tkmDisplayAreaRef this, 
               tkm_tOrientation  iOrientation );
DspA_tErr DspA_SetZoomLevel          ( tkmDisplayAreaRef this,
               int               inLevel );
DspA_tErr DspA_SetZoomCenter         ( tkmDisplayAreaRef this, 
               VoxelRef          ipCenter );
DspA_tErr DspA_SetZoomCenterToCursor ( tkmDisplayAreaRef this );
DspA_tErr DspA_HiliteSurfaceVertex   ( tkmDisplayAreaRef this,
               tkm_tSurfaceType  inSurface,
               int               inVertex );
DspA_tErr DspA_SetDisplayFlag        ( tkmDisplayAreaRef this,
               DspA_tDisplayFlag iWhichFlag,
               tBoolean          ibNewValue );
DspA_tErr DspA_ToggleDisplayFlag     ( tkmDisplayAreaRef this,
               DspA_tDisplayFlag iWhichFlag );
DspA_tErr DspA_SetTool               ( tkmDisplayAreaRef this,
               DspA_tTool        iTool );
DspA_tErr DspA_SetBrush              ( tkmDisplayAreaRef this,
               int               inRadius,
               DspA_tBrushShape  iShape,
               tBoolean          ib3D );
DspA_tErr DspA_SetBrushThreshold     ( tkmDisplayAreaRef this,
               tVolumeValue      inLow,
               tVolumeValue      inHigh,
               tVolumeValue      inNewValue );
DspA_tErr DspA_ChangeSliceBy_        ( tkmDisplayAreaRef this,
               int               inDelta );

/* only one display can be focused at a time. focusing on one will unfocus
   the previously focused one. */
DspA_tErr DspA_Focus    ( tkmDisplayAreaRef this );
DspA_tErr DspA_Blur     ( tkmDisplayAreaRef this );

/* get the size and location */
DspA_tErr DspA_GetPosition ( tkmDisplayAreaRef this,
           xPoint2nRef       opLocation,
           int*              onWidth,
           int*              onHeight );

/* routes events to specialized handlers */
DspA_tErr DspA_HandleEvent ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent );

/* internal handlers */
DspA_tErr DspA_HandleMouseUp_    ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent );
DspA_tErr DspA_HandleMouseDown_  ( tkmDisplayAreaRef this, 
          xGWin_tEventRef   ipEvent );
DspA_tErr DspA_HandleMouseMoved_ ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent );
DspA_tErr DspA_HandleKeyDown_    ( tkmDisplayAreaRef this, 
           xGWin_tEventRef   ipEvent );

/* uses the current brush settings to run the input function on a bunch
   of voxels based on the input point. */
DspA_tErr DspA_BrushVoxels_ ( tkmDisplayAreaRef this,
            VoxelRef          ipStartingVox,
            void(*ipFunction)(VoxelRef) );
void DspA_BrushVoxelsInThreshold_ ( VoxelRef ipVoxel );

/* for drawing the surface. */
tBoolean xUtil_FaceIntersectsPlane( MRI_SURFACE*     ipSurface,
            face_type*       ipFace,
            int              inPlane,
            tkm_tSurfaceType iSurface,
            tkm_tOrientation iOrientation );
void     xUtil_NormalizeVertexToVoxel( vertex_type*     ipVertex,
               tkm_tSurfaceType iSurface,
               tkm_tOrientation iOrientation,
               VoxelRef         opVoxel );
tBoolean xUtil_LineIntersectsPlane( VoxelRef         ipLineVoxA,
            VoxelRef         ipLineVoxB,
            int              inPlane,
            tBoolean         ibInterpolate,
            xPoint2fRef      opIntersectionPt );
DspA_tErr DspA_AdjustSurfaceDrawPoint_( tkmDisplayAreaRef this,
          xPoint2fRef       ipPoint );
DspA_tErr DspA_ParsePointList_( tkmDisplayAreaRef this,
        GLenum            inMode,
        xListRef          ipList );

/* schedule a redraw */
DspA_tErr DspA_Redraw_ ( tkmDisplayAreaRef this );

/* do the actual drawing */
DspA_tErr DspA_HandleDraw_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameBuffer_        ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSurface_            ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSurfaceDirect_      ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawCursor_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameAroundDisplay_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawAxes_               ( tkmDisplayAreaRef this );

/* build the frame buffer */
DspA_tErr DspA_BuildCurrentFrame_            ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFunctionalOverlayToFrame_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawControlPointsToFrame_     ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSelectionToFrame_         ( tkmDisplayAreaRef this );

/* other drawing subfunctions */
DspA_tErr DspA_BuildSurfaceDrawLists_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawCrosshairIntoFrame_ ( tkmDisplayAreaRef this,
           float*            ifaColor,
           xPoint2nRef       ipWhere,
           int               inSize );
void DspA_DrawVerticalArrow_ ( xPoint2nRef iStart,
             int         inLength,
             char*       isLabel );
void DspA_DrawHorizontalArrow_ ( xPoint2nRef iStart,
         int         inLength,
         char*       isLabel );

/* get info about the drawing state */
DspA_tErr DspA_GetCursor        ( tkmDisplayAreaRef this, 
          VoxelRef          opCursor );
DspA_tErr DspA_GetOrientation   ( tkmDisplayAreaRef this, 
          tkm_tOrientation* oOrientation );
DspA_tErr DspA_GetZoomLevel     ( tkmDisplayAreaRef this, 
          int*              oZoomLevel );
int DspA_GetCurrentSliceNumber_ ( tkmDisplayAreaRef this );

/* handles the lists of surface points */
DspA_tErr DspA_InitSurfaceLists_( tkmDisplayAreaRef this,
          int               inNumLists );
DspA_tErr DspA_PurgeSurfaceLists_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_NewSurfaceList_    ( tkmDisplayAreaRef this,
            tkm_tOrientation  iOrientation,
            tkm_tSurfaceType  iSurface,
            int               inSlice );
xListRef DspA_GetSurfaceList_       ( tkmDisplayAreaRef this,
              tkm_tOrientation  iOrientation,
              tkm_tSurfaceType  iSurface,
              int               inSlice );
int DspA_GetNumSurfaceLists_      ( tkmDisplayAreaRef this );
int DspA_GetSurfaceListIndex_     ( tkmDisplayAreaRef this,
            tkm_tOrientation  iOrientation,
            tkm_tSurfaceType  iSurface,
            int               inSlice );
DspA_tErr DspA_IsSurfaceCachced_( tkmDisplayAreaRef this,
          MRI_SURFACE*      ipSurface,
          tBoolean*         obIsCached );
DspA_tErr DspA_InitHashTables_( tkmDisplayAreaRef this,
        MRI_SURFACE*      ipSurface );
         

/* schedule a redraw */

/* conversion funcs */
DspA_tErr DspA_ConvertVolumeToBuffer_ ( tkmDisplayAreaRef this,
          VoxelRef          ipVolumeVox,
          xPoint2nRef       opBufferPt );

DspA_tErr DspA_ConvertBufferToVolume_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipBufferPt,
          VoxelRef          opVolumeVox );

DspA_tErr DspA_ConvertBufferToScreen_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipBufferPt,
          xPoint2nRef       opScreenPt );

DspA_tErr DspA_ConvertScreenToBuffer_ ( tkmDisplayAreaRef this,
          xPoint2nRef       opScreenPt,
          xPoint2nRef       ipBufferPt );

DspA_tErr DspA_ConvertPlaneToVolume_ ( tkmDisplayAreaRef this,
               xPoint2nRef       ipPlanePt,
               int               inSlice,
               tkm_tOrientation  iOrientation,
               VoxelRef          opVolumeVox );

/* send all viewing info to tcl */
DspA_tErr DspA_SendViewStateToTcl_ ( tkmDisplayAreaRef this );

DspA_tErr DspA_Verify       ( tkmDisplayAreaRef this );
DspA_tErr DspA_VerifyVolumeVoxel_  ( tkmDisplayAreaRef this,
             VoxelRef          ipVoxel );
DspA_tErr DspA_VerifyScreenPoint_ ( tkmDisplayAreaRef this,
             xPoint2nRef       ipScreenPt );
DspA_tErr DspA_VerifyBufferPoint_ ( tkmDisplayAreaRef this,
             xPoint2nRef       ipBufferPt );
/* set up opengl port */
void DspA_SetUpOpenGLPort_ ( tkmDisplayAreaRef this );


/* debug print */
void DspA_DebugPrint_ ( tkmDisplayAreaRef this );
void DspA_Signal ( char* isFuncName, int inLineNum, DspA_tErr ieCode );
char* DspA_GetErrorString ( DspA_tErr ieCode );


#endif








