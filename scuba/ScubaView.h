#ifndef ScubaView_h
#define ScubaView_h

#include "string_fixed.h"
#include <gl.h>
#include "View.h"
#include "DataCollection.h"
#include "InputState.h"
#include "ViewState.h"
#include "Layer.h"
#include "ScubaWindowToRASTranslator.h"
#include "ScubaToolState.h"
#include "ScubaTransform.h"
#include "Listener.h"

class ScubaView : public View, public ScubaWindowToRASTranslator, public Listener {

  friend class ScubaViewTester;

 public:

  ScubaView();
  virtual ~ScubaView();

  static int const kBytesPerPixel;         // for buffer size
  static int const kcInPlaneMarkerColors;  // number of preset colors

  // Sets the view. Used by something that wants to explicitly set up
  // the view area, such as a linked view broadcasting its position or
  // a script.
  void Set2DRASCenter ( float iRASCenter[3] );
  void Set2DZoomLevel ( float iZoom );
  void Set2DInPlane ( ViewState::Plane iPlane );

  void Get2DRASCenter ( float oRASCenter[3] );
  float Get2DZoomLevel ();
  ViewState::Plane Get2DInPlane ();

  // Add and remove layers that this view at a specific level. Note
  // that only one layer can be at a specific level, and adding a
  // layer at a level that is already occupied will replace the
  // existing layer.
  void SetLayerAtLevel ( int iLayerID, int iLevel );
  int GetLayerAtLevel ( int iLevel );
  void RemoveAllLayers ();
  void RemoveLayerAtLevel ( int iLevel );

  // Sets the same layers in another view.
  void CopyLayerSettingsToView ( ScubaView& iView );

  // Set the display transform for this view.
  void SetWorldToViewTransform ( int iTransformID );
  int GetWorldToViewTransform ();

  // Handle Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Handle broadcast messages.
  virtual void
    DoListenToMessage ( std::string isCommand, void* iData );

  // Implement ScubaWindowToRASTranslator.
  void TranslateWindowToRAS ( int iWindow[2], float oRAS[3] );
  void TranslateRASToWindow ( float iRAS[3], int oWindow[2] );

  // Get the first draw level with no layer assigned to it.
  int GetFirstUnusedDrawLevel ();

  // Set the flag to rebuild the draw overlay.
  void RebuildOverlayDrawList () { mbRebuildOverlayDrawList = true; }

  // Access the left/right flip flag.
  void SetFlipLeftRightYZ ( bool iFlip );
  bool GetFlipLeftRightYZ () { return mbFlipLeftRightInYZ; }

  // Get the inplane marker color.
  void GetInPlaneMarkerColor ( float oColor[3] );

protected:

  // Tells all the layers to draw in the correct order to the frame
  // buffer and then writes the frame buffer to the GL context.
  virtual void DoDraw ();

  // Resizes the frame buffer.
  virtual void DoReshape ( int iWidth, int iHeight );

  // Passes to layers.
  virtual void DoTimer ();
  
  // On mouse moves, this calls GetInfoAtRAS on all Layers and writes
  // the info as strings on the window.
  virtual void DoMouseMoved ( int iWindow[2], 
			      InputState& iInput, ScubaToolState& iTool );

  // Mouse up sets a marker at the current location. Mouse down and
  // mouse up may trigger tool effects.
  virtual void DoMouseUp ( int iWindow[2], 
			   InputState& iInput, ScubaToolState& iTool );
  virtual void DoMouseDown ( int iWindow[2], 
			     InputState& iInput, ScubaToolState& iTool );

  // Key up and down may trigger commands.
  virtual void DoKeyDown ( int iWindow[2], 
			   InputState& iInput, ScubaToolState& iTool );
  virtual void DoKeyUp ( int iWindow[2], 
			 InputState& Translates, ScubaToolState& iTool );

  // iInput window coords to RAS coordinates based on the current
  // view port.
  float ConvertWindowToRAS ( float iWindow, float iRASCenter, 
			     float iWindowDimension );
  float ConvertRASToWindow ( float iRAS, float iRASCenter, 
			     float iWindowDimension );

  // Sets a marker in the view.
  void SetMarker ();

  // The different steps in building our display. BuildFrameBuffer()
  // tells all the layers to copy their data to the frame
  // buffer. DrawFrameBuffer() copies it to the screen. BuildOverlay()
  // makes the draw lists for the gl command overlay. DrawOverlay()
  // executes the draw list.
  void BuildFrameBuffer ();
  void DrawFrameBuffer ();
  void BuildOverlay();
  void DrawOverlay ();

  void SetLinkedStatus ( bool ibLinked ) {mViewIDLinkedList[GetID()]=ibLinked;}
  bool GetLinkedStatus () { return mViewIDLinkedList[GetID()]; }

  // The draw list for the view overlay and a boolean saying whether
  // it should be rebuilt, usually when the view changes. This view
  // will actually use list kOverlayDrawListID + mID.
  #define kOverlayDrawListID 1
  bool mbRebuildOverlayDrawList;

  std::map<std::string,std::map<std::string,std::string> > mLabelValueMaps;

  // List of layers and their levels (level, layerID).
  std::map<int,int> mLevelLayerIDMap;
  
  // Current view information for this view.
  ViewState mViewState;

  // For keeping track of mouse movement navigation;
  float mMouseMoveDelta[2];
  int mLastMouseMoved[2];
  int mLastMouseDown[2];
  float mOriginalCenterRAS[3];
  float mOriginalZoom;

  float mInPlaneMovementIncrements[3];

  // The buffer for this view.
  GLubyte* mBuffer;

  // Key assignments.
  std::string msMoveViewLeft;
  std::string msMoveViewRight;
  std::string msMoveViewUp;
  std::string msMoveViewDown;
  std::string msMoveViewIn;
  std::string msMoveViewOut;
  std::string msZoomViewIn;
  std::string msZoomViewOut;

  // Link stuff.
  static int mCurrentBroadcaster;
  static std::map<int,bool> mViewIDLinkedList;

  // The world to view transform, or the 'view transform.' Applied to
  // all world coordinates before getting to the view.
  ScubaTransform* mWorldToView;

  // Whether to flip the right/left coordinates.
  bool mbFlipLeftRightInYZ;

  // The color to use when drawing this view's inplane on another view.
  float mInPlaneMarkerColor[3];

};  

class ScubaViewFactory : public ViewFactory {
 public:
  virtual View* NewView() { 
    return new ScubaView();
  }
};

class ScubaViewBroadcaster : public Broadcaster {
 public:
  static ScubaViewBroadcaster& GetBroadcaster ();
};



#endif
