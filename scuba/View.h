#ifndef View_h
#define View_h

#include <string>
#include "DebugReporter.h"

class View : public DebugReporter {

 public:
  View( std::string iID );
  ~View();

  void Draw();
  void Reshape();
  void Timer();
  void MouseMoved( int inX, int inY, int iButton, int iModifiers );
  void MouseUp( int inX, int inY, int iButton, int iModifers );
  void MouseDown( int inX, int inY, int iButton, int iModifers );
  void KeyDown( int inX, int inY, std::string isKey, int iModifers );
  void KeyUp( int inX, int inY, std::string isKey, int iModifers );

  void SetWidth( int iWidth ) { mWidth = iWidth; }
  void SetHeight( int iHeight ) { mHeight = iHeight; }

 protected:
  int mWidth;
  int mHeight;
  
  std::string msID;
};

#endif
