#ifndef tkmDisplayArea_h
#define tkmDisplayArea_h

#include "xTypes.h"
#include <GL/glut.h>
#include "xVoxel.h"
#include "tkmedit.h"
#include "xGLutWindow.h"
#include "xList.h"
#include "tkmFunctionalVolume.h"
#include "mriSurface.h"
#include "xGrowableArray.h"
#include "mriHeadPointList.h"
#include "mriROIGroup.h"

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
  DspA_tErr_ErrorAccessingSurface,
  DspA_tErr_OutOfMemory,
  DspA_tErr_ErrorAccessingControlPoints,
  DspA_tErr_ErrorAccessingSelection,
  DspA_tErr_ErrorAccessingWindow,
  DspA_tErr_ErrorAccessingFunctionalVolume,
  DspA_tErr_ErrorAccessingSurfaceList,
  DspA_tErr_ErrorAccessingHeadPointList,
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
  DspA_tDisplayFlag_ROIGroupOverlay,
  DspA_tDisplayFlag_FocusFrame,
  DspA_tDisplayFlag_Axes,
  DspA_tDisplayFlag_MaxIntProj,
  DspA_tDisplayFlag_HeadPoints,
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
  
  DspA_tMarker_Crosshair = 0,
  DspA_tMarker_Diamond,
  DspA_knNumMarkers
} DspA_tMarker;

#define DspA_kSignature 0x194ffb2

struct tkmDisplayArea {

  tSignature             mSignature;
  
  /* superpane info */
  struct tkmMeditWindow* mpWindow;
  
  /* our size and location */
  xPoint2n               mLocationInSuper;
  int                    mnWidth;
  int                    mnHeight;
  
  /* frame buffer */
  int                    mnVolumeSize;
  GLubyte*               mpFrameBuffer;
  float                  mfFrameBufferScaleX;
  float                  mfFrameBufferScaleY;

  /* view state */
  xVoxelRef              mpCursor;
  mri_tOrientation       mOrientation;
  int                    mnZoomLevel;
  xVoxelRef              mpZoomCenter;
  int                    mnHilitedVertexIndex;
  Surf_tVertexSet        mHilitedSurface;
  tBoolean               mabDisplayFlags [DspA_knNumDisplayFlags];
  tBoolean               mbSliceChanged;
  HPtL_tHeadPointRef     mpSelectedHeadPoint;
  int                    mnROIGroupIndex;
  
  /* surface lists */
  xGrowableArrayRef*     maSurfaceLists;

  /* display data */
  tVolumeRef             mpVolume;
  tVolumeRef             mpAuxVolume;
  mriROIGroupRef         mROIGroup;
  mriSurfaceRef          mpSurface;
  tkmFunctionalVolumeRef mpFunctionalVolume;
  x3DListRef             mpControlPoints;
  xListRef               mpSelectedControlPoints;
  x3DListRef             mpSelection;
  mriHeadPointListRef    mHeadPoints;
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
DspA_tErr DspA_SetROIGroup                   ( tkmDisplayAreaRef this,
                 mriROIGroupRef    iGroup );
DspA_tErr DspA_SetSurface                    ( tkmDisplayAreaRef this, 
                 mriSurfaceRef     ipSurface );
DspA_tErr DspA_SetOverlayVolume              ( tkmDisplayAreaRef this,
                 tkmFunctionalVolumeRef ipVol );
DspA_tErr DspA_SetControlPointsSpace         ( tkmDisplayAreaRef this,
                 x3DListRef        ipVoxels );
DspA_tErr DspA_SetControlPointsSelectionList ( tkmDisplayAreaRef this,
                 xListRef          ipVoxels );
DspA_tErr DspA_SetSelectionSpace             ( tkmDisplayAreaRef this, 
                 x3DListRef        ipVoxels );
DspA_tErr DspA_SetHeadPointList              ( tkmDisplayAreaRef this,
                 mriHeadPointListRef iList );

/* viewing state changes */
DspA_tErr DspA_SetCursor             ( tkmDisplayAreaRef this, 
              xVoxelRef          ipCursor );
DspA_tErr DspA_SetSlice              ( tkmDisplayAreaRef this,
               int               inSlice );
DspA_tErr DspA_SetOrientation        ( tkmDisplayAreaRef this, 
               mri_tOrientation  iOrientation );
DspA_tErr DspA_SetZoomLevel          ( tkmDisplayAreaRef this,
               int               inLevel );
DspA_tErr DspA_SetZoomCenter         ( tkmDisplayAreaRef this, 
              xVoxelRef          ipCenter );
DspA_tErr DspA_SetZoomCenterToCursor ( tkmDisplayAreaRef this );
DspA_tErr DspA_HiliteSurfaceVertex   ( tkmDisplayAreaRef this,
               Surf_tVertexSet  inSurface,
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
DspA_tErr DspA_ChangeSliceBy         ( tkmDisplayAreaRef this,
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
DspA_tErr DspA_GetSlice              ( tkmDisplayAreaRef this,
               int*              onSlice );

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
            xVoxelRef          ipStartingVox,
            void(*ipFunction)(xVoxelRef) );
void DspA_BrushVoxelsInThreshold_ (xVoxelRef ipVoxel );

/* select the currently clicked roi */
DspA_tErr DspA_SelectCurrentROI ( tkmDisplayAreaRef this );

/* for drawing the surface. */
tBoolean xUtil_FaceIntersectsPlane( MRI_SURFACE*     ipSurface,
            face_type*       ipFace,
            int              inPlane,
            Surf_tVertexSet iSurface,
            mri_tOrientation iOrientation );
void     xUtil_NormalizeVertexToVoxel( vertex_type*     ipVertex,
               Surf_tVertexSet iSurface,
               mri_tOrientation iOrientation,
              xVoxelRef         opVoxel );
tBoolean xUtil_LineIntersectsPlane(xVoxelRef         ipLineVoxA,
           xVoxelRef         ipLineVoxB,
            int              inPlane,
            tBoolean         ibInterpolate,
            xPoint2fRef      opIntersectionPt );
DspA_tErr DspA_AdjustSurfaceDrawPoint_( tkmDisplayAreaRef this,
          xPoint2fRef       ipPoint );
DspA_tErr DspA_ParsePointList_( tkmDisplayAreaRef this,
        GLenum            inMode,
        xGrowableArrayRef ipList );

/* schedule a redraw */
DspA_tErr DspA_Redraw_ ( tkmDisplayAreaRef this );

/* do the actual drawing */
DspA_tErr DspA_HandleDraw_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameBuffer_        ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSurface_            ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawHeadPoints_         ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawControlPoints_      ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawCursor_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameAroundDisplay_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawAxes_               ( tkmDisplayAreaRef this );

/* build the frame buffer */
DspA_tErr DspA_BuildCurrentFrame_            ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFunctionalOverlayToFrame_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSelectionToFrame_         ( tkmDisplayAreaRef this );

/* other drawing subfunctions */
DspA_tErr DspA_BuildSurfaceDrawLists_  ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawMarker_             ( tkmDisplayAreaRef this,
           DspA_tMarker      iType,
           float*            ifaColor,
           xPoint2nRef       ipWhere,
           int               inSize );
void DspA_DrawVerticalArrow_           ( xPoint2nRef iStart,
           int         inLength,
           char*       isLabel );
void DspA_DrawHorizontalArrow_         ( xPoint2nRef iStart,
           int         inLength,
           char*       isLabel );

/* get info about the drawing state */
DspA_tErr DspA_GetCursor        ( tkmDisplayAreaRef this, 
         xVoxelRef          opCursor );
DspA_tErr DspA_GetOrientation   ( tkmDisplayAreaRef this, 
          mri_tOrientation* oOrientation );
DspA_tErr DspA_GetZoomLevel     ( tkmDisplayAreaRef this, 
          int*              oZoomLevel );
int DspA_GetCurrentSliceNumber_ ( tkmDisplayAreaRef this );

/* tkmedit needs to get the selected head pt */
DspA_tErr DspA_GetSelectedHeadPt ( tkmDisplayAreaRef   this,
           HPtL_tHeadPointRef* opHeadPoint );

/* handles the lists of surface points */
DspA_tErr DspA_InitSurfaceLists_     ( tkmDisplayAreaRef this,
               int               inNumLists );
DspA_tErr DspA_PurgeSurfaceLists_    ( tkmDisplayAreaRef this );
DspA_tErr DspA_NewSurfaceList_       ( tkmDisplayAreaRef this,
               mri_tOrientation  iOrientation,
               Surf_tVertexSet   iSurface,
               int               inSlice );
xGrowableArrayRef DspA_GetSurfaceList_ ( tkmDisplayAreaRef this,
           mri_tOrientation  iOrientation,
           Surf_tVertexSet   iSurface,
           int               inSlice );
int DspA_GetNumSurfaceLists_         ( tkmDisplayAreaRef this );
int DspA_GetSurfaceListIndex_        ( tkmDisplayAreaRef this,
               mri_tOrientation  iOrientation,
               Surf_tVertexSet   iSurface,
               int               inSlice );
DspA_tErr DspA_IsSurfaceCachced_     ( tkmDisplayAreaRef this,
               MRI_SURFACE*      ipSurface,
               tBoolean*         obIsCached );

/* conversion funcs */
DspA_tErr DspA_ConvertVolumeToBuffer_ ( tkmDisplayAreaRef this,
          xVoxelRef         ipVolumeVox,
          xPoint2nRef       opBufferPt );

DspA_tErr DspA_ConvertBufferToVolume_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipBufferPt,
          xVoxelRef         opVolumeVox );

DspA_tErr DspA_ConvertBufferToScreen_ ( tkmDisplayAreaRef this,
          xPoint2nRef       ipBufferPt,
          xPoint2nRef       opScreenPt );

DspA_tErr DspA_ConvertScreenToBuffer_ ( tkmDisplayAreaRef this,
          xPoint2nRef       opScreenPt,
          xPoint2nRef       ipBufferPt );

DspA_tErr DspA_ConvertPlaneToVolume_ ( tkmDisplayAreaRef  this,
               xPoint2nRef        ipPlanePt,
               int                inSlice,
               mri_tOrientation   iOrientation,
               xVoxelRef          opVolumeVox );

/* send all viewing info to tcl */
DspA_tErr DspA_SendViewStateToTcl_ ( tkmDisplayAreaRef this );

DspA_tErr DspA_Verify       ( tkmDisplayAreaRef this );
DspA_tErr DspA_VerifyVolumexVoxl_  ( tkmDisplayAreaRef this,
            xVoxelRef          ipVoxel );
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








