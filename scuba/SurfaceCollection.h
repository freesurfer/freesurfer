#ifndef SurfaceCollection_H
#define SurfaceCollection_H

#include "string_fixed.h"
#include "DataCollection.h"
extern "C" {
#include "mrisurf.h"
#ifdef X
  #undef X 
#endif
#ifdef Y
  #undef Y
#endif
#ifdef Z
  #undef Z
#endif
}

class SurfaceCollection : public DataCollection {

 public:
  SurfaceCollection();
  virtual ~SurfaceCollection();

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "Surface"; }

  void SetSurfaceFileName ( std::string& ifnMRIS );

  void LoadSurface ();

  MRIS* GetMRIS();

  void LoadPatch ( std::string& ifnPatch );

  // For getting the surface to data transform from a volume.
  void SetDataToSurfaceTransformFromVolume ( VolumeCollection& iVolume );
  void SetDataToSurfaceTransformToDefault ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

  virtual ScubaROI* DoNewROI ();

  void RASToSurface  ( float iRAS[3], float oSurface[3] );
  void SurfaceToRAS  ( float iSurface[3], float oRAS[3] );

  // Surface access functions.
  int GetNumFaces ();
  int GetNumVerticesPerFace_Unsafe ( int inFace );
  void GetNthVertexInFace_Unsafe ( int inFace, int inVertex, 
				   float oRAS[3], bool* oRipped );

  int GetNumVertices ();
  void GetNthVertex_Unsafe ( int inVertex, 
			     float oRAS[3], bool* oRipped );

 protected:
  std::string mfnMRIS;
  MRIS* mMRIS;

  // We're dealing with two transforms here, data->world which is
  // defined in DataCollection, and surface->data which is defined
  // here. Actually it goes world->data->surface and
  // surface->data->world. We composite these in the last transform.
  Transform44 mDataToSurfaceTransform;
  Transform44 mWorldToSurfaceTransform;

  bool mbIsUsingVolumeForTransform;
  VolumeCollection* mTransformVolume;

  void CalcWorldToSurfaceTransform ();
};


#endif
