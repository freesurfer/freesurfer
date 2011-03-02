/**
 * @file  SurfaceCollection.h
 * @brief Implements DataCollection with a MRIS data set
 *
 * Very simple class that reads MRIS data and provides access to
 * geometry information.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.16 $
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


#ifndef SurfaceCollection_H
#define SurfaceCollection_H

#include "string_fixed.h"
#include "DataCollection.h"
extern "C" {
#include "mrisurf.h"
#include "mrishash.h"
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
  virtual std::string GetTypeDescription() {
    return "Surface";
  }

  void SetSurfaceFileName ( std::string& ifnMRIS );

  void LoadSurface ();

  MRIS* GetMRIS();

  void LoadPatch ( std::string& ifnPatch );

  // Return the bounds of the data in RAS coords. 0=xmin, 1=xmax,
  // 2=ymin, etc.
  virtual void GetDataRASBounds ( float oBounds[6] );

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

  int FindNearestVertexToRAS ( float iRAS[3], float *oMinDistance );

  int FindVertexAtRAS ( float iRAS[3], float *oMinDistance );

  // Surface access functions.
  int GetNumFaces ();
  int GetNumVerticesPerFace_Unsafe ( int inFace );
  void GetNthVertexInFace_Unsafe ( int inFace, int inVertex,
                                   float oRAS[3], bool* oRipped );

  int GetNumVertices ();
  void GetNthVertex_Unsafe ( int inVertex,
                             float oRAS[3], bool* oRipped );

  bool GetUseRealRAS ();

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

  // Hash table for finding vertices.
  MRIS_HASH_TABLE* mHashTable;

  // Bounds cache.
  bool mbBoundsCacheDirty;
  float mRASBounds[6];
};


#endif
