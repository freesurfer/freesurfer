#ifndef ScubaToolState_h
#define ScubaToolState_h

#include "IDTracker.h"
#include "TclCommandManager.h"

class ScubaToolState : public TclCommandListener, public IDTracker<ScubaToolState> {

 public:
  ScubaToolState();
  virtual ~ScubaToolState();

  enum Mode { navigation, plane, marker, voxelEditing, 
	      roiEditing, straightPath, edgePath };
  void SetMode ( Mode iMode ) { mMode = iMode; }
  Mode GetMode () { return mMode; }

  void SetLayerTarget ( int iTarget ) { mLayerTarget = iTarget; }
  int GetLayerTarget () { return mLayerTarget; }

  void SetNewValue ( float iValue ) { mNewValue = iValue; }
  float GetNewValue () { return mNewValue; }

  void SetOnlyBrushZero ( bool ibOnlyZero ) { mbOnlyFillZero = ibOnlyZero; }
  bool GetOnlyBrushZero () { return mbOnlyFillZero; }

  void SetBrushRadius ( float iRadius ) { mBrushRadius = iRadius; }
  float GetBrushRadius () { return mBrushRadius; }

  enum Shape { square, circle };
  void SetBrushShape ( Shape iShape ) { mBrushShape = iShape; }
  Shape GetBrushShape () { return mBrushShape; }

  void SetBrush3D ( bool ib3D ) { mbBrush3D = ib3D; }
  bool GetBrush3D () { return mbBrush3D; }

  void SetFloodStopAtPaths ( bool ibStop ) { mbFloodStopAtPaths = ibStop; }
  bool GetFloodStopAtPaths () { return mbFloodStopAtPaths; }

  void SetFloodStopAtROIs ( bool ibStop ) { mbFloodStopAtROIs = ibStop; }
  bool GetFloodStopAtROIs () { return mbFloodStopAtROIs; }

  void SetFloodFuzziness ( int iFuzziness ) { mFloodFuzziness = iFuzziness; }
  int GetFloodFuzziness () { return mFloodFuzziness; } 
  
  void SetFloodMaxDistance ( int iDistance ) { mFloodMaxDistance = iDistance; }
  int GetFloodMaxDistance () { return mFloodMaxDistance; } 
  
  void SetFlood3D ( bool ib3D ) { mbFlood3D = ib3D; }
  bool GetFlood3D () { return mbFlood3D; }

  void SetFloodSourceCollection ( int iCol ) { mFloodSourceCollection = iCol; }
  int GetFloodSourceCollection () { return mFloodSourceCollection; }

  void SetEdgePathStraightBias ( float iBias ) {mEdgePathStraightBias = iBias;}
  float GetEdgePathStraightBias () { return mEdgePathStraightBias; }

  void SetEdgePathEdgeBias ( float iBias ) {mEdgePathEdgeBias = iBias;}
  float GetEdgePathEdgeBias () { return mEdgePathEdgeBias; }

  void SetOnlyFloodZero ( bool ibOnlyZero ) { mbOnlyFloodZero = ibOnlyZero; }
  bool GetOnlyFloodZero () { return mbOnlyFloodZero; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

 protected:

  Mode mMode;

  // Layer target.
  int mLayerTarget;

  // Brush settings.
  float mBrushRadius;
  Shape mBrushShape;
  bool mbBrush3D;

  // Voxel edit settings.
  float mNewValue;
  float mbOnlyFillZero;

  // Flood settings.
  bool mbFloodStopAtPaths;
  bool mbFloodStopAtROIs;
  int mFloodFuzziness;
  int mFloodMaxDistance;
  bool mbFlood3D;
  int mFloodSourceCollection;
  bool mbOnlyFloodZero;

  // Edge path settings.
  float mEdgePathStraightBias;
  float mEdgePathEdgeBias;
};


#endif
