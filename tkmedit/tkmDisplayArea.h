#ifndef tkmDisplayArea_h
#define tkmDisplayArea_h

#include "xTypes.h"
#include <GL/glut.h>
#include "xVoxel.h"
#include "tkmedit.h"
#include "xGLutWindow.h"
#include "xList.h"
#include "mriVolume.h"
#include "tkmFunctionalVolume.h"
#include "mriSurface.h"
#include "mriColorLookupTable.h"
#include "xGrowableArray.h"
#include "mriHeadPointList.h"
#include "gca.h"
#include "const.h"
#include "vlabels.h"

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
  DspA_tErr_ErrorAccessingList,
  DspA_tErr_CouldntFindClosestVoxel,
  DspA_tErr_InvalidErrorCode,
  DspA_knNumErrorCodes
} DspA_tErr;

/* synced to tkm_interface.tcl */
typedef enum {
  DspA_tDisplayFlag_None = 0, 
  DspA_tDisplayFlag_AuxVolume,
  DspA_tDisplayFlag_Anatomical,
  DspA_tDisplayFlag_Cursor,
  DspA_tDisplayFlag_MainSurface,
  DspA_tDisplayFlag_OriginalSurface,
  DspA_tDisplayFlag_PialSurface,
  DspA_tDisplayFlag_InterpolateSurfaceVertices,
  DspA_tDisplayFlag_DisplaySurfaceVertices,
  DspA_tDisplayFlag_ControlPoints,
  DspA_tDisplayFlag_Selection,
  DspA_tDisplayFlag_FunctionalOverlay,
  DspA_tDisplayFlag_FunctionalColorScaleBar,
  DspA_tDisplayFlag_MaskToFunctionalOverlay,
  DspA_tDisplayFlag_HistogramPercentChange,
  DspA_tDisplayFlag_SegmentationVolumeOverlay,
  DspA_tDisplayFlag_AuxSegmentationVolume,
  DspA_tDisplayFlag_SegLabelVolumeCount,
  DspA_tDisplayFlag_DTIOverlay,
  DspA_tDisplayFlag_VectorField,
  DspA_tDisplayFlag_FocusFrame,
  DspA_tDisplayFlag_UndoVolume,
  DspA_tDisplayFlag_Axes,
  DspA_tDisplayFlag_MaxIntProj,
  DspA_tDisplayFlag_HeadPoints,
  DspA_knNumDisplayFlags
} DspA_tDisplayFlag;

/* synced to tkm_interface.tcl */
typedef enum {
  DspA_tTool_None         = -1,
  DspA_tTool_Navigate     =  0,
  DspA_tTool_SelectVoxels,
  DspA_tTool_EditVoxels,
  DspA_tTool_EditSegmentation,
  DspA_tTool_EditCtrlPts,
  DspA_knNumTools
} DspA_tTool;

/* synced to tkm_interface.tcl */
typedef enum {
  DspA_tBrushShape_None   = -1,
  DspA_tBrushShape_Circle =  0,
  DspA_tBrushShape_Square,
  DspA_knNumBrushShapes
} DspA_tBrushShape;

/* the number of different brush settings we have. they all share the
   same shape tho */
typedef enum {
  DspA_tBrush_None    = -1,
  DspA_tBrush_EditOne = 0,
  DspA_tBrush_EditTwo,
  DspA_knNumBrushes
} DspA_tBrush;

/* the threshold settings for a brush */
typedef struct {
  int mnLow;
  int mnHigh;
  int mnNewValue;
} DspA_tBrushInfo, *DspA_tBrushInfoRef;

typedef enum {
  DspA_tBrushTarget_None = -1,
  DspA_tBrushTarget_Main,
  DspA_tBrushTarget_MainAux,
  DspA_knNumBrushTargets
} DspA_tBrushTarget;

/* the combined brush settings. all use the same shape, but there can be
   different threshold settings (for different mouse buttons). */
typedef struct {
  
  /* target */
  DspA_tBrushTarget mTarget;

  /* shape */
  int mnRadius;
  DspA_tBrushShape mShape;
  tBoolean mb3D;
  
  /* thresholds */
  DspA_tBrushInfo mInfo[ DspA_knNumBrushes ];
  
} DspA_tBrushSettings;

/* segmentation brush settings. */
typedef struct {
  int               mNewValue;  /* The value that will be set */
  tBoolean          mb3D;
  tkm_tVolumeTarget mSrc;       /* Volume to use as source voxels. */
  int               mnFuzzy;
  int               mnDistance;
  tkm_tSegType      mDest;      /* The volume to affect; determined by which
				   seg volume is active. */

  int               mnPaintValue; /* Allows the brush to have two 'colors' */
  int               mnEraseValue; /* depending on which button is used. */

} DspA_tSegBrushSettings;


/* flood select settings. */
typedef struct {
  tBoolean           mb3D;
  tkm_tVolumeTarget  mSrc;       /* Volume to use as source voxels. */
  int                mnFuzzy;
  int                mnDistance;
} DspA_tFloodSelectSettings;

typedef enum {
  DspA_tSelectAction_None   = -1,
  DspA_tSelectAction_Select = 0,
  DspA_tSelectAction_Deselect,
  DspA_knNumSelectActions
} DspA_tSelectAction;

typedef enum {
  DspA_tMarker_None = -1,
  DspA_tMarker_Crosshair = 0,
  DspA_tMarker_Diamond,
  DspA_knNumMarkers
} DspA_tMarker;

typedef enum {
  DspA_tDisplaySet_Cursor = 0,
  DspA_tDisplaySet_Mouseover,
  DspA_knNumDisplaySets
} DspA_tDisplaySet;

/* parameters for the tcl histogram bar chart */
#define DspA_knHistoTitleLength 256
typedef struct {
  
  char    msTitle[DspA_knHistoTitleLength];
  char    msXAxisTitle[DspA_knHistoTitleLength];
  char    msYAxisTitle[DspA_knHistoTitleLength];
  char    msLabel1[DspA_knHistoTitleLength];
  char    msLabel2[DspA_knHistoTitleLength];
  int     mnNumValues;
  char**  masXAxisLabels;
  float*  mafValues1;
  float*  mafValues2;
  
} DspA_tHistogramParams, *DspA_tHistogramParamsRef;

/* entry in a surface list. has a boolean for whether it's a vertex or
   face, the original surface vertex (in ana idx), and the interpolated
   vertex. */
typedef struct {

  tBoolean mbVertex;

  /* the original vertex. */
  int      mnOriginalVertexIndex;
  xVoxel   mOriginalVertex;

  /* the neighboring vertex. */
  int      mnNeighborVertexIndex;
  xVoxel   mNeighborVertex;

  /* the original vertex projected onto the viewing plane. */
  xPoint2f mIntersectionPoint;
  
  /* the intersection of the edge between the original and neighboring
     vertices with the viewing plane. */
  xPoint2f mInterpIntersectionPoint;

  /* the interpIntersectionPoint in 3D space. */
  xVoxel   mInterpVertex;

  /* color for the vertex */
  tBoolean mOverrideColor;
  xColor3f mColor;

} DspA_tSurfaceListNode, *DspA_tSurfaceListNodeRef;

#define DspA_kSignature 0x194ffb2

struct tkmDisplayArea {
  
  tSignature             mSignature;
  
  /* superpane info */
  struct tkmMeditWindow* mpWindow;
  int                    mID;
  
  /* our size and location */
  xPoint2n               mLocationInSuper;
  int                    mnWidth;
  int                    mnHeight;
  
  /* frame buffer */
  int                    mnVolumeSizeX; /* This is actually the dimensions */
  int                    mnVolumeSizeY; /* of the 'screen space', which is */
  int                    mnVolumeSizeZ; /* 256^3. */
  GLubyte*               mpFrameBuffer;
  float                  mfFrameBufferScaleX;
  float                  mfFrameBufferScaleY;
  
  /* view state */
  xVoxelRef              mpLastCursor;
  xVoxelRef              mpCursor;
  xVoxelRef              mpMouseLocationAnaIdx;
  mri_tOrientation       mOrientation;
  int                    mnZoomLevel;
  xVoxelRef              mpZoomCenter;
  int                    mnHilitedVertexIndex;
  Surf_tVertexSet        mHilitedSurface;
  tBoolean               mabDisplayFlags [DspA_knNumDisplayFlags];
  tBoolean               mbSliceChanged;
  HPtL_tHeadPointRef     mpSelectedHeadPoint;
  int                    mnSegmentationVolumeIndex;
  int        manSurfaceLineWidth[tkm_knNumSurfaceTypes][Surf_knNumVertexSets];
  xColor3f   maSurfaceLineColor[tkm_knNumSurfaceTypes][Surf_knNumVertexSets];

  float                  mfSegmentationAlpha;
  tkm_tAxis              maDTIAxisForComponent[xColr_knNumComponents];
  float                  mfDTIAlpha;

  /* for navigation tool */
  xVoxelRef              mpOriginalZoomCenter;
  int                    mnOriginalSlice;
  int                    mnOriginalZoomLevel;
  xPoint2n               mLastClick;
  xPoint2f               mTotalDelta;
  
  /* updated whenever we get a mouse moved event */
  xPoint2n               mMouseLocation;
  
  /* surface lists */
  xGrowableArrayRef*     maSurfaceLists[tkm_knNumSurfaceTypes];
  
  /* display data */
  mriVolumeRef           mpVolume[tkm_knNumVolumeTypes];
  mriVolumeRef           mSegmentationVolume[tkm_knNumSegTypes];
  mriColorLookupTableRef mSegmentationColorTable[tkm_knNumSegTypes];
  mriSurfaceRef          mpSurface[tkm_knNumSurfaceTypes];
  tkmFunctionalVolumeRef mpFunctionalVolume;
  x3DListRef             mpControlPoints;
  mriVolumeRef           mpSelection;
  mriHeadPointListRef    mHeadPoints;
  GCA*                   mGCAVolume;
  TRANSFORM*             mGCATransform;
  VLI*                   mVLI1 ;
  VLI*                   mVLI2 ;
  char                   isVLI1_name[STRLEN] ;
  char                   isVLI2_name[STRLEN] ;
  mriVolumeRef           mpDTIVolume;
};
typedef struct tkmDisplayArea tkmDisplayArea;
typedef tkmDisplayArea *tkmDisplayAreaRef;

#define DspA_knNumBytesPerPixel     4
#define DspA_knMaxPixelValue        (GLubyte)255
#define DspA_knRedPixelCompIndex    0
#define DspA_knGreenPixelCompIndex  1
#define DspA_knBluePixelCompIndex   2
#define DspA_knAlphaPixelCompIndex  3

#define DspA_knMaxZoomLevel 128
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

/* Super will set this */
DspA_tErr DspA_SetID ( tkmDisplayAreaRef this, int inID );

/* update window title */
DspA_tErr DspA_UpdateWindowTitle ( tkmDisplayAreaRef this );

/* all the things it will draw */
DspA_tErr DspA_SetVolume                     ( tkmDisplayAreaRef this,
                                               mriVolumeRef      ipVolume,
                                               int               inSizeX,
                                               int               inSizeY,
                                               int               inSizeZ );
DspA_tErr DspA_SetAuxVolume                  ( tkmDisplayAreaRef this,
                                               mriVolumeRef      ipVolume,
                                               int               inSizeX,
                                               int               inSizeY,
                                               int               inSizeZ) ;
DspA_tErr DspA_SetSegmentationVolume         ( tkmDisplayAreaRef this,
					       tkm_tSegType      iType,
					       mriVolumeRef      iVolume );
DspA_tErr DspA_SetSegmentationColorTable     ( tkmDisplayAreaRef this,
					       tkm_tSegType      iType,
					       mriColorLookupTableRef iCLUT );
DspA_tErr DspA_SetSurface                    ( tkmDisplayAreaRef this, 
					       tkm_tSurfaceType  iType,
					       mriSurfaceRef     ipSurface );
DspA_tErr DspA_SetOverlayVolume              ( tkmDisplayAreaRef this,
					       tkmFunctionalVolumeRef ipVol );
DspA_tErr DspA_SetControlPointsSpace         ( tkmDisplayAreaRef this,
					       x3DListRef        ipVoxels );
DspA_tErr DspA_SetSelectionSpace             ( tkmDisplayAreaRef this, 
					       mriVolumeRef      ipVolume );
DspA_tErr DspA_SetHeadPointList              ( tkmDisplayAreaRef this,
					       mriHeadPointListRef iList );
DspA_tErr DspA_SetGCA                        ( tkmDisplayAreaRef this,
					       GCA*              iVolume,
					       TRANSFORM*        iTransform );
DspA_tErr DspA_SetVLIs                       ( tkmDisplayAreaRef this,
					       VLI*              iVLI1,
					       VLI*              iVLI2,
					       char*             isVLI1_name,
					       char*             isVLI2_name);
DspA_tErr DspA_SetDTIVolume                  ( tkmDisplayAreaRef  this, 
					       mriVolumeRef       iVolume );

/* viewing state changes */
DspA_tErr DspA_SetCursor             ( tkmDisplayAreaRef this, 
				       xVoxelRef         ipCursor );
DspA_tErr DspA_ConvertAndSetCursor   ( tkmDisplayAreaRef this,
				       mri_tCoordSpace   iFromSpace,
				       xVoxelRef         ipCoord );
DspA_tErr DspA_SetSlice              ( tkmDisplayAreaRef this,
				       int               inSlice );
DspA_tErr DspA_SetOrientation        ( tkmDisplayAreaRef this, 
				       mri_tOrientation  iOrientation );
DspA_tErr DspA_SetZoomLevel          ( tkmDisplayAreaRef this,
				       int               inLevel );
DspA_tErr DspA_SetZoomCenter         ( tkmDisplayAreaRef this, 
				       xVoxelRef         ipCenter );
DspA_tErr DspA_SetZoomCenterToCursor ( tkmDisplayAreaRef this );
DspA_tErr DspA_HiliteSurfaceVertex   ( tkmDisplayAreaRef this,
				       Surf_tVertexSet   inSurface,
				       int               inVertex );
DspA_tErr DspA_SetDisplayFlag        ( tkmDisplayAreaRef this,
				       DspA_tDisplayFlag iWhichFlag,
				       tBoolean          ibNewValue );
DspA_tErr DspA_ToggleDisplayFlag     ( tkmDisplayAreaRef this,
				       DspA_tDisplayFlag iWhichFlag );
DspA_tErr DspA_SetTool               ( tkmDisplayAreaRef this,
				       DspA_tTool        iTool );
DspA_tErr DspA_SetBrushTarget        ( tkmDisplayAreaRef this,
				       DspA_tBrushTarget iTarget );
DspA_tErr DspA_SetBrushShape         ( tkmDisplayAreaRef this,
				       int               inRadius,
				       DspA_tBrushShape  iShape,
				       tBoolean          ib3D );
DspA_tErr DspA_SetBrushInfo          ( tkmDisplayAreaRef this,
				       DspA_tBrush       iBrush,
				       DspA_tBrushInfoRef iInfo ); 
DspA_tErr DspA_SetCursorColor        ( tkmDisplayAreaRef this,
				       xColor3fRef       iColor );
DspA_tErr DspA_SetCursorShape        ( tkmDisplayAreaRef this,
				       DspA_tMarker      iShape );
DspA_tErr DspA_SetSurfaceLineWidth   ( tkmDisplayAreaRef this,
				       tkm_tSurfaceType  iType,
				       Surf_tVertexSet   iSurface,
				       int               inWidth );
DspA_tErr DspA_SetSurfaceLineColor   ( tkmDisplayAreaRef this,
				       tkm_tSurfaceType  iType,
				       Surf_tVertexSet   iSurface,
				       xColor3fRef       iColor );
DspA_tErr DspA_SetFloodSelectParams ( tkmDisplayAreaRef this,
				       DspA_tFloodSelectSettings* iSettings );
DspA_tErr DspA_SetSegBrushInfo      ( tkmDisplayAreaRef this,
				       DspA_tSegBrushSettings* iSettings );
DspA_tErr DspA_ChangeSliceBy         ( tkmDisplayAreaRef this,
				       int               inDelta );
DspA_tErr DspA_SetSegmentationAlpha  ( tkmDisplayAreaRef this,
				       float             ifAlpha );
DspA_tErr DspA_SetDTIAlpha           ( tkmDisplayAreaRef this,
				       float             ifAlpha );
DspA_tErr DspA_SetDTIAxisForComponent ( tkmDisplayAreaRef this,
					tkm_tAxis         iAxis,
					xColr_tComponent  iComponent );


/* only one display can be focused at a time. focusing on one will unfocus
   the previously focused one. */
DspA_tErr DspA_Focus    ( tkmDisplayAreaRef this );
DspA_tErr DspA_Blur     ( tkmDisplayAreaRef this );

/* get the size and location */
DspA_tErr DspA_GetPosition ( tkmDisplayAreaRef this,
			     xPoint2nRef       opLocation,
			     int*              onWidth,
			     int*              onHeight );
DspA_tErr DspA_GetSlice    ( tkmDisplayAreaRef this,
			     int*              onSlice );

/* Sets the cursor to the center of the selection volume. */
DspA_tErr DspA_SetCursorToCenterOfSpace ( tkmDisplayAreaRef this,
					  mriVolumeRef      ipVolume );

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

/* set brush to default values */
DspA_tErr DspA_SetBrushInfoToDefault ( tkmDisplayAreaRef this,
				       DspA_tBrush       iBrush );

/* uses the current brush settings to run the input function on a
   bunch of voxels based on the input point. The callback function
   takes an array of xVoxelRef*, a length of the array, and a pointer
   to data. */
DspA_tErr DspA_BrushVoxels_ ( tkmDisplayAreaRef this,
			      xVoxelRef         ipStartingVox,
			      void*             ipData,
			      void(*ipFunction)(xVoxelRef,int,void*) );
/* parameter here is a DspA_tBrush */
void DspA_BrushVoxelsInThreshold_ ( xVoxelRef ipaVoxel, int inCount, void* );
/* parameter here is a DspA_tSelectAction */
void DspA_SelectVoxels_           ( xVoxelRef ipaVoxel, int inCount, void* );
/* no parameter here */
void DspA_EditSegmentationVoxels_ ( xVoxelRef ipaVoxel, int inCount, void* );

/* select the currently clicked seg label */
DspA_tErr DspA_SelectCurrentSegLabel ( tkmDisplayAreaRef this );

/* graph the avg of the currently clicked seg label */
DspA_tErr DspA_GraphCurrentSegLabelAvg ( tkmDisplayAreaRef this );

/* gets current mouse position, translates to volume idx, and sends
   info to tcl as mouseover info. */
DspA_tErr DspA_SendMouseInfoToTcl ( tkmDisplayAreaRef this );

/* for drawing the surface. */
/* these functions convert a vertex in any orientation to a normalized
   orientation for comparison and then back again. used before and
   after a voxel goes into the xUtil_LineIntersectsPlane function. */
DspA_tErr DspA_NormalizeVoxel_   ( xVoxelRef        ipAnaIdx,
				   mri_tOrientation iOrientation,
				   xVoxelRef        opNormIdx );
DspA_tErr DspA_UnnormalizeVoxel_ ( xVoxelRef        ipNormIdx,
				   mri_tOrientation iOrientation,
				   xVoxelRef        opAnaIdx );

/* this sees if two voxels intersect a plane in 3D and if so, returns
   true, as well as the 2D intersection points on that plane; one the
   simple projection of the first voxel onto the plane, and the other
   the actual point on the line between the two voxels where it
   interesects the plane. */
tBoolean xUtil_LineIntersectsPlane ( xVoxelRef         ipAnaIdxA,
				     xVoxelRef         ipAnaIdxB,
				     int               inPlane,
				     xPoint2fRef       opIntersectionPt,
				     xPoint2fRef     opInterpIntersectionPt );

/* walks through a surafce list and draws the points with openGL. can
   be used to draw lines or points. the port should already be set up
   with the proper color and point size. */
DspA_tErr DspA_ParsePointList_ ( tkmDisplayAreaRef this,
				 GLenum            inMode,
				 xGrowableArrayRef ipList );

/* used to move a drawn crosshair point from .0,.0,.0 to .5,.5,.5 to
   make it show in the middle of a voxel. */
DspA_tErr DspA_AdjustSurfaceAnaIdx   ( tkmDisplayAreaRef this,
				       xVoxelRef         iAnaIdx );
DspA_tErr DspA_UnadjustSurfaceAnaIdx ( tkmDisplayAreaRef this,
				       xVoxelRef         iAnaIdx );

/* schedule a redraw */
DspA_tErr DspA_Redraw_ ( tkmDisplayAreaRef this );

/* do the actual drawing */
DspA_tErr DspA_HandleDraw_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameBuffer_        ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSurface_            ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawHeadPoints_         ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawControlPoints_      ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawVectorField_        ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawCursor_             ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFrameAroundDisplay_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawAxes_               ( tkmDisplayAreaRef this );

/* build the frame buffer */
DspA_tErr DspA_BuildCurrentFrame_              ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSegmentationOverlayToFrame_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawDTIOverlayToFrame_          ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawUndoableVoxelsOverlayToFrame_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawFunctionalOverlayToFrame_   ( tkmDisplayAreaRef this );
DspA_tErr DspA_DrawSelectionToFrame_           ( tkmDisplayAreaRef this );

/* other drawing subfunctions */
DspA_tErr DspA_BuildSurfaceDrawLists_  ( tkmDisplayAreaRef this,
					 tkm_tSurfaceType  iType );
DspA_tErr DspA_DrawMarker_             ( tkmDisplayAreaRef this,
					 DspA_tMarker      iType,
					 float*            ifaColor,
					 xPoint2nRef       ipWhere,
					 int               inSize );
DspA_tErr DspA_DrawVector_             ( tkmDisplayAreaRef this,
					 float*            ifaColor,
					 xVoxelRef         ipVoxelStart,
					 xVoxelRef         ipVoxelDirection );
void      DspA_DrawVerticalArrow_      ( xPoint2nRef       iStart,
					 int               inLength,
					 char*             isLabel );
void      DspA_DrawHorizontalArrow_    ( xPoint2nRef       iStart,
					 int               inLength,
					 char*             isLabel );

/* get info about the drawing state */
DspA_tErr DspA_GetCursor        ( tkmDisplayAreaRef this, 
				  xVoxelRef          opCursor );
DspA_tErr DspA_GetOrientation   ( tkmDisplayAreaRef this, 
				  mri_tOrientation* oOrientation );
DspA_tErr DspA_GetZoomLevel     ( tkmDisplayAreaRef this, 
				  int*              oZoomLevel );
int DspA_GetCurrentSliceNumber_ ( tkmDisplayAreaRef this );

/* draw data into the histogram window */
DspA_tErr DspA_DrawHistogram ( tkmDisplayAreaRef        this,
			       DspA_tHistogramParamsRef iParams );

/* write the current distance to the vertex closest to the cursor */
DspA_tErr DspA_SetSurfaceDistanceAtCursor ( tkmDisplayAreaRef this );

/* looks at which volume edge the cursor is nearest and then 'cuts' the
   region from that edge to the same plane at the cursor by setting
   those values to 0.*/
DspA_tErr DspA_SmartCutAtCursor ( tkmDisplayAreaRef this );

/* tkmedit needs to get the selected head pt */
DspA_tErr DspA_GetSelectedHeadPt ( tkmDisplayAreaRef   this,
				   HPtL_tHeadPointRef* opHeadPoint );

/* find the closest interpolated surface point. */
DspA_tErr DspA_GetClosestInterpSurfVoxel ( tkmDisplayAreaRef this,
					   tkm_tSurfaceType  iType,
					   Surf_tVertexSet   iSet,
					   xVoxelRef         iAnaIdx,
					   xVoxelRef         oOrigAnaIdx,
					   xVoxelRef         oInterpAnaIdx,
					   char*             osDescription );

/* handles the lists of surface points */
DspA_tErr DspA_InitSurfaceLists_     ( tkmDisplayAreaRef this,
				       int               inNumLists );
DspA_tErr DspA_PurgeSurfaceLists_    ( tkmDisplayAreaRef this );
DspA_tErr DspA_NewSurfaceList_       ( tkmDisplayAreaRef this,
				       tkm_tSurfaceType  iType,
				       mri_tOrientation  iOrientation,
				       Surf_tVertexSet   iSurface,
				       int               inSlice );
xGrowableArrayRef DspA_GetSurfaceList_ ( tkmDisplayAreaRef this,
					 tkm_tSurfaceType  iType,
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
DspA_tErr DspA_ConvertVolumeToBufferf_ ( tkmDisplayAreaRef this,
					 xVoxelRef         ipVolumeVox,
					 xPoint2fRef       opBufferPt );

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
				       xPoint2fRef        ipPlanePt,
				       int                inSlice,
				       mri_tOrientation   iOrientation,
				       xVoxelRef          opVolumeVox );

DspA_tErr DspA_ConvertVolumeToPlane_ ( tkmDisplayAreaRef this,
				       xVoxelRef         ipVolumeVox,
				       mri_tOrientation  iOrientation,
				       xPoint2fRef       opPlanePt,
				       int*              onSlice );

/* send all viewing info to tcl */
DspA_tErr DspA_SendViewStateToTcl_ ( tkmDisplayAreaRef this );
DspA_tErr DspA_SendPointInformationToTcl_ ( tkmDisplayAreaRef this,
					    DspA_tDisplaySet  iSet,
					    xVoxelRef         iAnaIdx );

DspA_tErr DspA_Verify       ( tkmDisplayAreaRef this );
DspA_tErr DspA_VerifyVolumeVoxel_ ( tkmDisplayAreaRef this,
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
