#ifndef ScubaLayer2DMRIS_h
#define ScubaLayer2DMRIS_h

#include "Layer.h"
#include "SurfaceCollection.h"

class ScubaLayer2DMRIS : public Layer {

  friend class ScubaLayer2DMRISTester;

 public:
  ScubaLayer2DMRIS ();
  virtual ~ScubaLayer2DMRIS ();

  // Associate a surface collection with this layer.
  void SetSurfaceCollection ( SurfaceCollection& iSurface );

  // Tell the layer to draw its contents into a GL frame buffer.
  virtual void DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				ViewState& iViewState,
				ScubaWindowToRASTranslator& iTranslator );
  
  // Asks the layer to describe a point of data by adding pairs of
  // labels and values.
  virtual void GetInfoAtRAS ( float inX, float inY, float inZ,
			   std::map<std::string,std::string>& iLabelValues );
  
  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () { return "2DMRIS"; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

 protected:
  SurfaceCollection* mSurface;

};


#endif
