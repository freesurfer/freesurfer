/**
 * @file  ScubaLayer2DMRIS.h
 * @brief Draws a MRIS surface
 *
 * Draws the intersection of a surface and a plane into the
 * Layer. Also handles getting information about a surface at an RAS
 * point, and responding to Tcl commands requesting vertex
 * information.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.18 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


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

  // Tell the layer to draw its contents that need openGl commands.
  virtual void DrawIntoGL ( ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator );

  // Asks the layer to describe a point of data by making InfoAtRAS
  // structs.
  virtual void GetInfoAtRAS ( float iRAS[3],
                              std::list<InfoAtRAS>& ioInfo );

  // These are the names of the reportable info values we provide.
  enum ReportableInfo { Vertex = 0, Distance, kcReportableInfo };
  static char* const kaReportableInfo[kcReportableInfo];

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription () {
    return "2DMRIS";
  }

  virtual void DataChanged();

  // Handle tools.
  virtual void HandleTool ( float iRAS[3], ViewState& iViewState,
                            ScubaWindowToRASTranslator& iTranslator,
                            ScubaToolState& iTool, InputState& iInput );

  void FindRASLocationOfVertex ( int inVertex, float oRAS[3] );

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

  // Return the surface collection.
  virtual DataCollection* GetMainDataCollection() {
    return mSurface;
  }

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

  void SetLineWidth ( int iWidth ) {
    mLineWidth = iWidth;
  }
  int GetLineWidth () {
    return mLineWidth;
  }

  void SetDrawVertices ( bool ibDrawVertices ) {
    mbDrawVertices = ibDrawVertices;
  }
  bool GetDrawVertices () {
    return mbDrawVertices;
  }

  // To process command line options.
  void ProcessOption ( std::string isOption, std::string isValue );

protected:
  SurfaceCollection* mSurface;

  int maLineColor[3];
  int maVertexColor[3];
  int mLineWidth;

  bool mbDrawVertices;

  std::list<int> mCachedDrawList;
  ViewState mCachedViewState;

  void ClearCache();
};


#endif
