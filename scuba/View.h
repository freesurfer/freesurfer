#ifndef View_h
#define View_h

#include <string>
#include "IDTracker.h"
#include "DebugReporter.h"

// This is a base class used by ScubaFrame as an interface to views,
// or the internal components of the Frame. Each Frame can have
// multiple non-overlapping Views, and each View can have multiple
// layers.

class View : public DebugReporter, public IDTracker<View> {
  
  friend class ViewTester;

public:
  View();
  virtual ~View();

  void Draw();
  void Reshape( int iWidth, int iHeight );
  void Timer();
  void MouseMoved( int inX, int inY, int iButton, int iModifiers );
  void MouseUp( int inX, int inY, int iButton, int iModifers );
  void MouseDown( int inX, int inY, int iButton, int iModifers );
  void KeyDown( int inX, int inY, std::string isKey, int iModifers );
  void KeyUp( int inX, int inY, std::string isKey, int iModifers );

  void SetWidth( int iWidth ) { mWidth = iWidth; }
  void SetHeight( int iHeight ) { mHeight = iHeight; }

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

 protected:
  virtual void DoDraw();
  virtual void DoReshape( int iWidth, int iHeight );
  virtual void DoTimer();
  virtual void DoMouseMoved( int inX, int inY, int iButton, int iModifiers );
  virtual void DoMouseUp( int inX, int inY, int iButton, int iModifers );
  virtual void DoMouseDown( int inX, int inY, int iButton, int iModifers );
  virtual void DoKeyDown( int inX, int inY, std::string isKey, int iModifers );
  virtual void DoKeyUp( int inX, int inY, std::string isKey, int iModifers );

  int mWidth;
  int mHeight;
  
  std::string msLabel;
};


class ViewFactory {
 public:
  virtual View* NewView() { 
    return new View();
  }
};



#endif
