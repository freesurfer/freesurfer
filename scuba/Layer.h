#ifndef Layer_h
#define Layer_h

#include <string>
#include <gl.h>
#include <map>
#include "ViewState.h"
#include "IDTracker.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "ScubaWindowToRASTranslator.h"

class LayerStaticTclListener : public DebugReporter, public TclCommandListener {
  public :
    ~LayerStaticTclListener ();

    virtual TclCommandResult
      DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );
};


class Layer : public DebugReporter, public IDTracker<Layer>, public TclCommandListener {

  friend class ScubaViewTester;

 public:

  Layer();
  virtual ~Layer();

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer( GLubyte* iBuffer, int iWidth, int iHeight,
			       ViewState& iViewState,
			       ScubaWindowToRASTranslator& iTranslator );
  
  // Asks the layer to describe a point of data by adding pairs of
  // labels and values.
  virtual void GetInfoAtRAS ( float inX, float inY, float inZ,
			   std::map<std::string,std::string>& iLabelValues );
  
  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "BaseLayer"; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

  void SetOpacity( float iOpacity ) { mOpacity = iOpacity; }
  float GetOpacity() { return mOpacity; }

  void SetWidth( int iWidth ) { mWidth = iWidth; }
  void SetHeight( int iHeight ) { mHeight = iHeight; }

 protected:

  int mWidth;
  int mHeight;

  std::string msLabel;
  
  float mOpacity;

  static bool mbRegisteredStaticListener;
  static LayerStaticTclListener mStaticListener;
};



#endif
