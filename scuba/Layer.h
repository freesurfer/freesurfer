#ifndef Layer_h
#define Layer_h

#include <string>
#include <gl.h>
#include <map>
#include "ViewState.h"
#include "IDTracker.h"
#include "DebugReporter.h"

class Layer : public DebugReporter, public IDTracker<Layer> {

 public:

  Layer();
  virtual ~Layer();

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer( GLbyte* iBuffer, ViewState& iViewState );
  
  // Asks the layer to describe a point of data by adding pairs of
  // labels and values.
  virtual void GetInfoAtRAS ( float inX, float inY, float inZ,
			   std::map<std::string,std::string>& iLabelValues );
  
  void SetLabel( std::string isLabel ) { msLabel = isLabel; }
  std::string GetLabel() { return msLabel; }

 protected:

  std::string msLabel;
  
  float mOpacity;
};



#endif
