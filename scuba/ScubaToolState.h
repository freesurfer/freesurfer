#ifndef ScubaToolState_h
#define ScubaToolState_h

#include "IDTracker.h"
#include "TclCommandManager.h"

class ScubaToolState : public TclCommandListener, public IDTracker<ScubaToolState> {

 public:
  ScubaToolState();
  virtual ~ScubaToolState();

  enum Mode { navigation, plane, marker, voxelEditing, 
	      roiEditing, straightLine, edgeLine };
  void SetMode ( Mode iMode ) { mMode = iMode; }
  Mode GetMode () { return mMode; }

  void SetBrushRadius ( float iRadius ) { mBrushRadius = iRadius; }
  float GetBrushRadius () { return mBrushRadius; }

  enum Shape { square, circle };
  void SetBrushShape ( Shape iShape ) { mBrushShape = iShape; }
  Shape GetBrushShape () { return mBrushShape; }

  void SetBrush3D ( bool ib3D ) { mbBrush3D = ib3D; }
  bool GetBrush3D () { return mbBrush3D; }

  void SetFloodStopAtLines ( bool ibStop ) { mbFloodStopAtLines = ibStop; }
  bool GetFloodStopAtLines () { return mbFloodStopAtLines; }

  void SetFloodStopAtROIs ( bool ibStop ) { mbFloodStopAtROIs = ibStop; }
  bool GetFloodStopAtROIs () { return mbFloodStopAtROIs; }

  void SetFloodFuzziness ( int iFuzziness ) { mFloodFuzziness = iFuzziness; }
  int GetFloodFuzziness () { return mFloodFuzziness; } 
  
  void SetFloodMaxDistance ( int iDistance ) { mFloodMaxDistance = iDistance; }
  int GetFloodMaxDistance () { return mFloodMaxDistance; } 
  
  void SetFlood3D ( bool ib3D ) { mbFlood3D = ib3D; }
  bool GetFlood3D () { return mbFlood3D; }

  void SetEdgeLineStraightBias ( float iBias ) {mEdgeLineStraightBias = iBias;}
  float GetEdgeLineStraightBias () { return mEdgeLineStraightBias; }

  void SetEdgeLineEdgeBias ( float iBias ) {mEdgeLineEdgeBias = iBias;}
  float GetEdgeLineEdgeBias () { return mEdgeLineEdgeBias; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

 protected:

  Mode mMode;

  // Brush settings.
  float mBrushRadius;
  Shape mBrushShape;
  bool mbBrush3D;

  // Flood settings.
  bool mbFloodStopAtLines;
  bool mbFloodStopAtROIs;
  int mFloodFuzziness;
  int mFloodMaxDistance;
  bool mbFlood3D;

  // Edge line settings.
  float mEdgeLineStraightBias;
  float mEdgeLineEdgeBias;
};


#endif
