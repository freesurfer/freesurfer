#ifndef ScubaToolState_h
#define ScubaToolState_h

#include "IDTracker.h"
#include "TclCommandManager.h"

class ScubaToolState : public TclCommandListener, public IDTracker<ScubaToolState> {

 public:
  ScubaToolState();
  virtual ~ScubaToolState();

  enum Mode { navigation, plane, marker, voxelEditing, voxelFilling,
	      roiEditing, roiFilling, straightPath, edgePath };
  enum Shape { voxel, square, circle };
  enum FuzzinessType { seed, gradient };

  // Accessors.
  void SetMode ( Mode iMode ) { mMode = iMode; }
  Mode GetMode () { return mMode; }

  // Target.
  void SetTargetLayer ( int iTarget ) { mTargetLayer = iTarget; }
  int GetTargetLayer () { return mTargetLayer; }

  // voxelEditing
  void SetNewValue ( float iValue ) { mNewValue = iValue; }
  float GetNewValue () { return mNewValue; }

  void SetEraseValue ( float iValue ) { mEraseValue = iValue; }
  float GetEraseValue () { return mEraseValue; }

  // Brush settings.
  void SetOnlyBrushZero ( bool ibOnlyZero ) { mbOnlyFillZero = ibOnlyZero; }
  bool GetOnlyBrushZero () { return mbOnlyFillZero; }

  void SetBrushRadius ( float iRadius ) { mBrushRadius = iRadius; }
  float GetBrushRadius () { return mBrushRadius; }

  void SetBrushShape ( Shape iShape ) { mBrushShape = iShape; }
  Shape GetBrushShape () { return mBrushShape; }
  std::string GetBrushShapeAsString ();

  void SetBrush3D ( bool ib3D ) { mbBrush3D = ib3D; }
  bool GetBrush3D () { return mbBrush3D; }

  // Fill settings
  void SetFloodStopAtPaths ( bool ibStop ) { mbFloodStopAtPaths = ibStop; }
  bool GetFloodStopAtPaths () { return mbFloodStopAtPaths; }

  void SetFloodStopAtROIs ( bool ibStop ) { mbFloodStopAtROIs = ibStop; }
  bool GetFloodStopAtROIs () { return mbFloodStopAtROIs; }

  void SetFloodFuzziness ( float iFuzziness ) { mFloodFuzziness = iFuzziness; }
  float GetFloodFuzziness () { return mFloodFuzziness; } 
  
  void SetFloodMaxDistance ( float iDistance ) { mFloodMaxDistance =iDistance;}
  float GetFloodMaxDistance () { return mFloodMaxDistance; } 
  
  void SetFlood3D ( bool ib3D ) { mbFlood3D = ib3D; }
  bool GetFlood3D () { return mbFlood3D; }

  void SetFloodSourceCollection ( int iCol ) { mFloodSourceCollection = iCol; }
  int GetFloodSourceCollection () { return mFloodSourceCollection; }

  void SetOnlyFloodZero ( bool ibOnlyZero ) { mbOnlyFloodZero = ibOnlyZero; }
  bool GetOnlyFloodZero () { return mbOnlyFloodZero; }

  void SetFuzzinessType ( FuzzinessType iType ) { mFloodFuzzinessType =iType; }
  FuzzinessType GetFuzzinessType () { return mFloodFuzzinessType; }

  // edgePath.
  void SetEdgePathStraightBias ( float iBias ) {mEdgePathStraightBias = iBias;}
  float GetEdgePathStraightBias () { return mEdgePathStraightBias; }

  void SetEdgePathEdgeBias ( float iBias ) {mEdgePathEdgeBias = iBias;}
  float GetEdgePathEdgeBias () { return mEdgePathEdgeBias; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

 protected:

  Mode mMode;

  // Layer target.
  int mTargetLayer;

  // Brush settings.
  float mBrushRadius;
  Shape mBrushShape;
  bool mbBrush3D;

  // Voxel edit settings.
  float mNewValue;
  float mEraseValue;
  float mbOnlyFillZero;

  // Flood settings.
  bool mbFloodStopAtPaths;
  bool mbFloodStopAtROIs;
  float mFloodFuzziness;
  float mFloodMaxDistance;
  bool mbFlood3D;
  int mFloodSourceCollection;
  bool mbOnlyFloodZero;
  FuzzinessType mFloodFuzzinessType;

  // Edge path settings.
  float mEdgePathStraightBias;
  float mEdgePathEdgeBias;
};


#endif
