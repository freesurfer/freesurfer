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
  virtual void GetInfoAtRAS ( float iRAS[3],
			   std::map<std::string,std::string>& iLabelValues );
  
  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () { return "2DMRIS"; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  void SetLineColor3d ( int iaLineColor[3] ) {
    maLineColor[0] = iaLineColor[0];
    maLineColor[1] = iaLineColor[1];
    maLineColor[2] = iaLineColor[2];
  }
  void GetLineColor3d ( int oaLineColor[3] ) {
    oaLineColor[0] = maLineColor[0];
    oaLineColor[1] = maLineColor[1];
    oaLineColor[2] = maLineColor[2];
  }
  void SetVertexColor3d ( int iaVertexColor[3] ) {
    maVertexColor[0] = iaVertexColor[0];
    maVertexColor[1] = iaVertexColor[1];
    maVertexColor[2] = iaVertexColor[2];
  }
  void GetVertexColor3d ( int oaVertexColor[3] ) {
    oaVertexColor[0] = maVertexColor[0];
    oaVertexColor[1] = maVertexColor[1];
    oaVertexColor[2] = maVertexColor[2];
  }


 protected:
  SurfaceCollection* mSurface;

  int maLineColor[3];
  int maVertexColor[3];
};


#endif
