#ifndef View_h
#define View_h

#include <string>
#include "FrameView.h"
#include "DataCollection.h"

class ViewLayer {

 public:
  typedef int ID;

  ViewLayer( ID );
  virtual ~ViewLayer();

  // Tell the layer to draw its contents into a GL frame buffer.
  void DrawIntoView( GLbyte* iBuffer, View::Plane iInPlane,  
		     float iXBegin, float iXEnd, float iXSpacing,
		     float iYBegin, float iYEnd, float iTSpacing,
		     float iZBegin, float iZEnd, float iZSpacing );

  // Asks the layer to describe a point of data by adding pairs of
  // labels and values.
  void GetInfoAtRAS ( float inX, float inY, float inZ,
		      std::map<std::string,std::string> iLabelValues );


  // Access to the layers by ID and to the whole ID list. 
  static ViewLayer& GetLayer( ID const iID );
  static std::list<ID>& GetLayerIDList ();

 protected:
  virtual void DoDraw();
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int inZ, 
			     int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int inZ,
			  int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int inZ,
			    int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, int inZ,
			  std::string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, int inZ,
			std::string isKey, int iModifers );

  float mOpacity;
  
  // All the allocated layers.
  static list<ID> mLayers;
};


class View : public FrameView {

 public:
  enum Plane { X, Y, Z };

  View( std::string iID );
  virtual ~View();

  // Sets the view. Used by something that wants to explicitly set up
  // the view area, such as a linked view broadcasting its position or
  // a script.
  void Set2DRASCenter ( float iRASCenter[] );
  void Set2DZoomLevel ( float iZoom );

  void Set3DEye ();
  void Set3DLook ();
  void Set3DUp ();

  // Sets a marker in the view.
  void SetMarker ();

  // Set the layers that this view renders.
  void SetLayers ();

protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, std::string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, std::string isKey, int iModifers );

  // Translates window coords to RAS coordinates based on the current
  // view port.
  void TranslateWindowToRAS ( int iXWindow, int iYWindow,
			      int& oXRAS, int& oYRAS, int& oZRAS );

  // List of layers.
  std::list<ViewLayer::ID> mLayers;
  
  // Orientation information for 2D views.
  Plane mInPlane;

  // Link information.
  std::list<View*> mLinkedViews;
};  

#endif
