/**
 * @file  vtkFSSurfaceSource.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/05/11 18:41:19 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef __vtkFSSurfaceSource_h
#define __vtkFSSurfaceSource_h

#include <vector>

#include "vtkSource.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
extern "C" {
#include "mrisurf.h"
}

class vtkFSSurfaceSource : public vtkSource {
public:

  static vtkFSSurfaceSource *New();
  vtkTypeRevisionMacro(vtkFSSurfaceSource,vtkSource);

  // Description:
  // This will call the MRISread function from the FS library.
  void MRISRead ( char const* ifnSurface );

  // Description:
  // Get the output of this source.
  vtkPolyData* GetOutput ();
  vtkPolyData* GetOutput ( int inOutput );
  void         SetOutput ( vtkPolyData* iOutput );

  // Description:
  // Coordinate conversion. RAS space is defined by various header
  // metadata and used by the Layer to display things in the right
  // space.
  void ConvertSurfaceToRAS ( float iSurfX, float iSurfY, float iSurfZ,
                             float& oRASX, float& oRASY, float& oRASZ ) const;
  void ConvertRASToSurface ( float iRASX, float iRASY, float iRASZ,
                             float& oSurfX, float& oSurfY, float& oSurfZ) const;
  void ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const;
  void ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const;

  void GetRASBounds ( float ioBounds[6] );

  vtkTransform const* GetSurfaceToRASTransform () const;

  float GetRASCenterX () const;
  float GetRASCenterY () const;
  float GetRASCenterZ () const;

  // Description:
  // Returns the number of vertices in the surface.
  int GetNumberOfVertices () const;

  // Description:
  // Get the vertex number from a RAS or surface RAS point. This uses
  // the hash table and finds only the closest vertex point. If
  // oDistance is not NULL, the distance to the found point will be
  // returned there.
  int FindVertexAtRAS        ( float const iRAS[3],        float* oDistance );
  int FindVertexAtSurfaceRAS ( float const iSurfaceRAS[3], float* oDistance );

  // Description:
  // Get the RAS or surface RAS coords at a vertex index.
  void GetRASAtVertex        ( int inVertex, float ioRAS[3] );
  void GetSurfaceRASAtVertex ( int inVertex, float ioRAS[3] );

  // Description:
  // Return a list a of vertex indices that form the shortest path
  // from inStartVertex to inEndVertex.
  void FindPath ( int inStartVertex, int inEndVertex, 
		  std::vector<int>& iolPath );

  MRIS* GetMRIS() { return mMRIS; }

protected:

  vtkFSSurfaceSource();
  ~vtkFSSurfaceSource();

  void Execute();

  // The surface.
  MRIS* mMRIS;

  //   [  0  1  2  3 ]
  //   [  4  5  6  7 ]
  //   [  8  9 10 11 ]
  //   [ 12 13 14 15 ]
  double mSurfaceToRASMatrix[16];
  vtkTransform* mSurfaceToRASTransform;
  
  // RAS bounds.
  bool mbBoundsCacheDirty;
  float mRASBounds[6];

  // The center in RAS space.
  float mRASCenter[3];

  // Hash table so we can look up vertices. Uses v->x,y,z.
  MRIS_HASH_TABLE* mHashTable;

private:
  vtkFSSurfaceSource(const vtkFSSurfaceSource&);  // Not implemented.
  void operator=(const vtkFSSurfaceSource&);  // Not implemented.
};

#endif
