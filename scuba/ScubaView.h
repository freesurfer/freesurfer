#ifndef ScubaView_h
#define ScubaView_h

#include <string>
#include <gl.h>
#include "View.h"
#include "DataCollection.h"
#include "InputState.h"
#include "ViewState.h"
#include "Layer.h"

class ScubaView : public View {

  friend class ScubaViewTester;

 public:

  ScubaView();
  virtual ~ScubaView();

  static int const kBytesPerPixel;

  // Sets the view. Used by something that wants to explicitly set up
  // the view area, such as a linked view broadcasting its position or
  // a script.
  void Set2DRASCenter ( float iRASCenter[] );
  void Set2DZoomLevel ( float iZoom );
  void Set2DInPlane ( ViewState::Plane iPlane );

  // Adds a layers that this view at a specific level. Note that only
  // one layer can be at a specific level, and adding a layer at a
  // level that is already occupied will replace the existing layer.
  void AddLayer ( int iLayerID, int iLevel );

  virtual void DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

protected:

  // Tells all the layers to draw in the correct order to the frame
  // buffer and then writes the frame buffer to the GL context.
  virtual void DoDraw();

  // Resizes the frame buffer.
  virtual void DoReshape( int iWidth, int iHeight );

  // Passes to layers.
  virtual void DoTimer();
  
  // On mouse moves, this calls GetInfoAtRAS on all Layers and writes
  // the info as strings on the window.
  virtual void DoMouseMoved( int inX, int inY, InputState& iState );

  // Mouse up sets a marker at the current location. Mouse down and
  // mouse up may trigger tool effects.
  virtual void DoMouseUp( int inX, int inY, InputState& iState );
  virtual void DoMouseDown( int inX, int inY, InputState& iState );

  // Key up and down may trigger commands.
  virtual void DoKeyDown( int inX, int inY, InputState& iState );
  virtual void DoKeyUp( int inX, int inY, InputState& Translates );

  // iState window coords to RAS coordinates based on the current
  // view port.
  void TranslateWindowToRAS ( int iXWindow, int iYWindow,
			      float& oXRAS, float& oYRAS, float& oZRAS );
  float ConvertWindowToRAS ( float iWindow, float iRASCenter, 
			     float iWindowDimension );

  // Sets a marker in the view.
  void SetMarker ();

  // List of layers and their levels (level, layerID).
  std::map<int,int> mLevelLayerIDMap;
  
  // Current view information for this view.
  ViewState mViewState;

  // The buffer for this view.
  GLbyte* mBuffer;
};  

class ScubaViewFactory : public ViewFactory {
 public:
  virtual View* NewView() { 
    return new ScubaView();
  }
};





#endif
