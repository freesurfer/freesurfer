#ifndef tkmMeditWindow_h
#define tkmMeditWindow_h

#include "tcl.h"
#include "xTypes.h"
#include "xDebug.h"
#include "xGLutWindow.h"
#include "tkmedit.h"
#include "tkmDisplayArea.h"
#include "tkmFunctionalVolume.h"
#include "mriSurface.h"
#include "mriHeadPointList.h"
#include "mriROIGroup.h"

typedef enum {

  MWin_tErr_NoErr = 0,
  MWin_tErr_AllocationFailed,
  MWin_tErr_DisplayAreaAllocationFailed,
  MWin_tErr_ErrorAccessingWindow,
  MWin_tErr_InvalidPtr,
  MWin_tErr_InvalidSignature,
  MWin_tErr_InvalidDisplayIndex,
  MWin_tErr_InvalidDisplayArea,
  MWin_tErr_InvalidDisplayConfiguration,
  MWin_tErr_InvalidCoordinateSpace,
  MWin_tErr_ErrorAccessingDisplay,
  MWin_tErr_WrongNumberArgs,
  MWin_tErr_InvalidEvent,
  MWin_tErr_InvalidErrorCode,
  MWin_knNumErrorCodes

} MWin_tErr;

typedef enum {

  MWin_tLinkPolicy_None = 0,
  MWin_tLinkPolicy_MultipleOrientations,
  MWin_tLinkPolicy_Mosaic,
  MWin_knNumLinkPolicies
} MWin_tLinkPolicy;

#define MWin_kSignature 0xffb21001

#define MWin_knMaxNumAreas     16
#define MWin_kAllDisplayAreas -1

#define MWin_knSpaceBetweenWindowAndPanel 20

struct tkmMeditWindow {

  tSignature       mSignature;

  /* window data */
  xGLutWindowRef   mpGLutWindow;
  int              mnID;
  int              mnWidth;
  int              mnHeight;

  /* our display areas */
  int                    mnRows;
  int                    mnCols;
  MWin_tLinkPolicy    mLinkPolicy;
  struct tkmDisplayArea* mapDisplays [MWin_knMaxNumAreas];
  int                    mnLastClickedArea;

  /* linking the cursor? */
  tBoolean         mbLinkedCursor;

  /* accepting commands from tcl? */
  tBoolean         mbAcceptingTclCommands;
};
typedef struct tkmMeditWindow tkmMeditWindow;
typedef tkmMeditWindow *tkmMeditWindowRef;

MWin_tErr MWin_New    ( tkmMeditWindowRef* oppWindow,
      char*              isTitle,
      int                inWidth, 
      int                inHeight );
MWin_tErr MWin_Delete ( tkmMeditWindowRef* ioppWindow );


/* set window title */
MWin_tErr MWin_SetWindowTitle ( tkmMeditWindowRef ipWindow,
        char*             isTitle );

/* get window width and height */
MWin_tErr MWin_GetWindowSize ( tkmMeditWindowRef ipWindow,
             int*              onX,
             int*              onY,
             int*              onWidth,
             int*              onHeight );

/* configuration the display set up */
MWin_tErr MWin_SetDisplayConfiguration ( tkmMeditWindowRef   this,
           int                 inCols,
           int                 inRows,
           MWin_tLinkPolicy iPolicy );
MWin_tErr MWin_PositionDisplays_       ( tkmMeditWindowRef this );

/* setting display info. specify by display area index, starting from 0
   for the upper-left area. use -1 to specify all areas. */
MWin_tErr MWin_SetVolume                     ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 tVolumeRef        ipVolume,
                 int               inSize );
MWin_tErr MWin_SetAuxVolume                  ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 tVolumeRef        ipVolume,
                 int               inSize );
MWin_tErr MWin_SetROIGroup                   ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 mriROIGroupRef    iGroup );
MWin_tErr MWin_SetSurface                    ( tkmMeditWindowRef this, 
                 int               inDispIndex,
                 mriSurfaceRef     ipSurface );
MWin_tErr MWin_SetOverlayVolume              ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 tkmFunctionalVolumeRef ipVol );
MWin_tErr MWin_SetControlPointsSpace         ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 x3DListRef        ipVoxels );
MWin_tErr MWin_SetControlPointsSelectionList ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 xListRef          ipVoxels );
MWin_tErr MWin_SetSelectionSpace             ( tkmMeditWindowRef this, 
                 int               inDispIndex,
                 x3DListRef        ipVoxels );
MWin_tErr MWin_SetHeadPointList              ( tkmMeditWindowRef this,
                 int               inDispIndex,
                 mriHeadPointListRef iList );


/* viewing state changes. specify the display area the same way as above. */
MWin_tErr MWin_SetLinkedCursorFlag   ( tkmMeditWindowRef this, 
               tBoolean          ibLinkCursor );
MWin_tErr MWin_ToggleLinkedCursorFlag( tkmMeditWindowRef this );
MWin_tErr MWin_SetCursor             ( tkmMeditWindowRef this, 
               int               inDispIndex,
               xVoxelRef         ipCursor );
MWin_tErr MWin_SetOrientation        ( tkmMeditWindowRef this, 
               int               inDispIndex,
               mri_tOrientation  iOrientation );
MWin_tErr MWin_SetZoomCenter         ( tkmMeditWindowRef this, 
               int               inDispIndex,
               xVoxelRef         ipCenter );
MWin_tErr MWin_SetZoomLevel          ( tkmMeditWindowRef this, 
               int               inDispIndex,
               int               inLevel );
MWin_tErr MWin_SetZoomCenterToCursor ( tkmMeditWindowRef this,
               int               inDispIndex );
MWin_tErr MWin_HiliteSurfaceVertex   ( tkmMeditWindowRef this,
               int               inDispIndex,
               Surf_tVertexSet  inSurface,
               int               inVertex );
MWin_tErr MWin_SetDisplayFlag        ( tkmMeditWindowRef this,
               int               inDispIndex,
               DspA_tDisplayFlag iWhichFlag,
               tBoolean          ibNewValue );
MWin_tErr MWin_SetVolumeColorScale   ( tkmMeditWindowRef this,
               int               inDispIndex,
               int               inMin,
               int               inMid,
               int               inMax );

/* get the viewing state of the last clicked display area */
MWin_tErr MWin_GetCursor         ( tkmMeditWindowRef   this,
           xVoxelRef           opCursor );
MWin_tErr MWin_GetOrientation    ( tkmMeditWindowRef   this,
           mri_tOrientation*   oOrientation );

/* tkmedit needs to get the selected head pt. bad design. */
MWin_tErr MWin_GetSelectedHeadPt ( tkmMeditWindowRef   this,
           HPtL_tHeadPointRef* opHeadPoint );

/* for cursor linking. a display area whose cursor was set calls this
   function. if we have cursor linking turned on, this will set all display
   cursors or flags. */
MWin_tErr MWin_CursorChanged       ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             xVoxelRef         ipCursor );
MWin_tErr MWin_ZoomLevelChanged    ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             int               inZoomLevel );
MWin_tErr MWin_DisplayFlagChanged  ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             DspA_tDisplayFlag iWhichFlag,
             tBoolean          ibNewValue );
MWin_tErr MWin_OrientationChanged  ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             mri_tOrientation  iOrientation );
MWin_tErr MWin_SliceChanged        ( tkmMeditWindowRef this,
             tkmDisplayAreaRef ipDisplay,
             int               inDelta );

/* callback for xGLutWindow. this just passes the ptr, a tkmMeditWindowRef,
   and the event to our normal event handler below. */
void MWin_EventCallback ( void*           ipWindow,
        xGWin_tEventRef ipEvent );

/* routes events to the display areas */
void MWin_HandleEvent ( tkmMeditWindowRef this, 
      xGWin_tEventRef   ipEvent );


/* takes redraw requests from tkmedit */
MWin_tErr MWin_Redraw ( tkmMeditWindowRef this );

/* this not only schedules a redraw but forces all display areas to 
   rebuild their slices. */
MWin_tErr MWin_RedrawAll ( tkmMeditWindowRef this );

/* do the actual drawing */
MWin_tErr MWin_HandleDraw_ ( tkmMeditWindowRef this );

/* move the tool window directly under the medit window */
MWin_tErr MWin_PlaceToolWindow_ ( tkmMeditWindowRef this );

/* focus on a different display area */
MWin_tErr MWin_ChangeFocusedDisplayAreaBy_( tkmMeditWindowRef this,
              int               inDelta );  

/* register tcl commands */
MWin_tErr MWin_RegisterTclCommands ( tkmMeditWindowRef this,
             Tcl_Interp*       ipInterp );

/* accept or reject tcl commands. this is done to keep the inital var
   settings in the tk window from affecting the medit window. */
MWin_tErr MWin_AcceptTclCommands ( tkmMeditWindowRef this,
           tBoolean          ibAccept );

MWin_tErr MWin_Verify             ( tkmMeditWindowRef this );
MWin_tErr MWin_VerifyDisplayIndex ( tkmMeditWindowRef this,
            int               nDisplayIndex );

char* MWin_GetErrorString ( MWin_tErr ieCode );

/* these tcl commands call the cooresponding on the display area that
   was last clicked on. */
int MWin_TclSetLinkedCursorFlag ( ClientData  iClientData, 
          Tcl_Interp* ipInterp,
          int         argc,
          char*       argv[] );
int MWin_TclSetCursor        ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetOrientation   ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetZoomLevel     ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetZoomCenter    ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetDisplayConfig ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetDisplayFlag   ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetTool          ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetBrush         ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclSetBrushThreshold ( ClientData  iClientData, 
        Tcl_Interp* ipInterp,
        int         argc,
        char*       argv[] );
int MWin_TclSelectCurrentROI ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );
int MWin_TclRedrawAll        ( ClientData  iClientData, 
             Tcl_Interp* ipInterp,
             int         argc,
             char*       argv[] );

#endif
