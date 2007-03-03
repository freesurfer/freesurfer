/**
 * @file  vtkFSSurfaceSource.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/03 00:04:11 $
 *    $Revision: 1.2 $
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


#include <stdexcept>
#include "vtkFSSurfaceSource.h"
#include "vtkObjectFactory.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"

using namespace std;

vtkStandardNewMacro( vtkFSSurfaceSource );
vtkCxxRevisionMacro( vtkFSSurfaceSource, "$Revision: 1.2 $" );

vtkFSSurfaceSource::vtkFSSurfaceSource() :
    mMRIS( NULL ),
    mSurfaceToRASTransform( NULL ),
    mbBoundsCacheDirty( true ),
    mHashTable( NULL ) {

  this->vtkSource::SetNthOutput(0, vtkPolyData::New());

  // Releasing data for pipeline parallism.
  // Filters will know it is empty.
  this->Outputs[0]->ReleaseData();
  this->Outputs[0]->Delete();
}

vtkFSSurfaceSource::~vtkFSSurfaceSource () {

  if ( NULL != mHashTable )
    MHTfree( &mHashTable );

}

void
vtkFSSurfaceSource::MRISRead( char const* ifn ) {

  char* fn = strdup( ifn );
  mMRIS = ::MRISread( fn );
  if ( mMRIS == NULL ) {
    throw runtime_error( "MRISread failed" );
  }

  // Get some info from the MRIS.
  if ( mMRIS->vg.valid ) {

    mRASCenter[0] = mMRIS->vg.c_r;
    mRASCenter[1] = mMRIS->vg.c_a;
    mRASCenter[2] = mMRIS->vg.c_s;

    mSurfaceToRASMatrix[0] = mMRIS->vg.x_r;
    mSurfaceToRASMatrix[1] = mMRIS->vg.y_r;
    mSurfaceToRASMatrix[2] = mMRIS->vg.z_r;
    mSurfaceToRASMatrix[3] = mMRIS->vg.c_r;
    mSurfaceToRASMatrix[4] = mMRIS->vg.x_a;
    mSurfaceToRASMatrix[5] = mMRIS->vg.y_a;
    mSurfaceToRASMatrix[6] = mMRIS->vg.z_a;
    mSurfaceToRASMatrix[7] = mMRIS->vg.c_a;
    mSurfaceToRASMatrix[8] = mMRIS->vg.x_s;
    mSurfaceToRASMatrix[9] = mMRIS->vg.y_s;
    mSurfaceToRASMatrix[10] = mMRIS->vg.z_s;
    mSurfaceToRASMatrix[11] = mMRIS->vg.c_s;
    mSurfaceToRASMatrix[12] = 0;
    mSurfaceToRASMatrix[13] = 0;
    mSurfaceToRASMatrix[14] = 0;
    mSurfaceToRASMatrix[15] = 1;
    
  } else if ( mMRIS->lta ) {

    mRASCenter[0] = mMRIS->lta->xforms[0].src.c_r;
    mRASCenter[1] = mMRIS->lta->xforms[0].src.c_a;
    mRASCenter[2] = mMRIS->lta->xforms[0].src.c_s;

    mSurfaceToRASMatrix[0] = 1;
    mSurfaceToRASMatrix[1] = 0;
    mSurfaceToRASMatrix[2] = 0;
    mSurfaceToRASMatrix[3] = -mRASCenter[0];
    mSurfaceToRASMatrix[4] = 0;
    mSurfaceToRASMatrix[5] = 1;
    mSurfaceToRASMatrix[6] = 0;
    mSurfaceToRASMatrix[7] = -mRASCenter[1];
    mSurfaceToRASMatrix[8] = 0;
    mSurfaceToRASMatrix[9] = 0;
    mSurfaceToRASMatrix[10] = 1;
    mSurfaceToRASMatrix[11] = -mRASCenter[2];
    mSurfaceToRASMatrix[12] = 0;
    mSurfaceToRASMatrix[13] = 0;
    mSurfaceToRASMatrix[14] = 0;
    mSurfaceToRASMatrix[15] = 1;

  }

  // Make our transform object and set the matrix.
  mSurfaceToRASTransform = vtkTransform::New();
  mSurfaceToRASTransform->SetMatrix( mSurfaceToRASMatrix );
  
  // Make the hash table. This makes it with v->x,y,z.
  if ( NULL != mHashTable )
    MHTfree( &mHashTable );
  mHashTable = MHTfillVertexTableRes( mMRIS, NULL, CURRENT_VERTICES, 2.0 );
}

vtkPolyData*
vtkFSSurfaceSource::GetOutput () {

  if ( this->NumberOfOutputs < 1 ) {
    return NULL;
  }

  return (vtkPolyData *)(this->Outputs[0]);
}

vtkPolyData*
vtkFSSurfaceSource::GetOutput ( int inOutput ) {

  return (vtkPolyData *)this->vtkSource::GetOutput( inOutput );
}

void
vtkFSSurfaceSource::SetOutput ( vtkPolyData* iOutput ) {

  this->vtkSource::SetNthOutput( 0, iOutput );
}

void
vtkFSSurfaceSource::ConvertSurfaceToRAS ( float iX, float iY, float iZ,
					  float& oX, float& oY, float& oZ ) const {

  float surface[3];
  float ras[3];

  surface[0] = iX;
  surface[1] = iY;
  surface[2] = iZ;

  this->ConvertSurfaceToRAS( surface, ras );

  oX = ras[0];
  oY = ras[1];
  oZ = ras[2];
}

void
vtkFSSurfaceSource::ConvertRASToSurface ( float iX, float iY, float iZ,
					  float& oX, float& oY, float& oZ ) const {

  float ras[3];
  float surface[3];
  
  ras[0] = iX;
  ras[1] = iY;
  ras[2] = iZ;
  
  this->ConvertRASToSurface( ras, surface );

  oX = surface[0];
  oY = surface[1];
  oZ = surface[2];
}

void
vtkFSSurfaceSource::ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const {

  mSurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void
vtkFSSurfaceSource::ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const {

  mSurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void
vtkFSSurfaceSource::GetRASBounds ( float oRASBounds[6] ) {

  if ( NULL == mMRIS ) {

    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
      oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;

    return;
  }

  if ( mbBoundsCacheDirty ) {

    mRASBounds[0] = mRASBounds[2] = mRASBounds[4] = 999999;
    mRASBounds[1] = mRASBounds[3] = mRASBounds[5] = -999999;

    // Find the bounds.
    for ( int vno = 0; vno < mMRIS->nvertices; vno++ ) {

      // Translate to actual RAS coords.
      float rasX, rasY, rasZ;
      this->ConvertSurfaceToRAS( mMRIS->vertices[vno].x,
                                 mMRIS->vertices[vno].y,
                                 mMRIS->vertices[vno].z,
                                 rasX, rasY, rasZ );

      if ( rasX < mRASBounds[0] ) mRASBounds[0] = rasX;
      if ( rasX > mRASBounds[1] ) mRASBounds[1] = rasX;
      if ( rasY < mRASBounds[2] ) mRASBounds[2] = rasY;
      if ( rasY > mRASBounds[3] ) mRASBounds[3] = rasY;
      if ( rasZ < mRASBounds[4] ) mRASBounds[4] = rasZ;
      if ( rasZ > mRASBounds[5] ) mRASBounds[5] = rasZ;

    }

    mbBoundsCacheDirty = false;
  }

  oRASBounds[0] = mRASBounds[0];
  oRASBounds[1] = mRASBounds[1];
  oRASBounds[2] = mRASBounds[2];
  oRASBounds[3] = mRASBounds[3];
  oRASBounds[4] = mRASBounds[4];
  oRASBounds[5] = mRASBounds[5];
}

vtkTransform const*
vtkFSSurfaceSource::GetSurfaceToRASTransform () const {

  return mSurfaceToRASTransform;
}

float
vtkFSSurfaceSource::GetRASCenterX () const {

  return mRASCenter[0];
}

float
vtkFSSurfaceSource::GetRASCenterY () const {

  return mRASCenter[1];
}

float
vtkFSSurfaceSource::GetRASCenterZ () const {

  return mRASCenter[2];
}

int
vtkFSSurfaceSource::GetNumberOfVertices () const {

  if( mMRIS )
    return mMRIS->nvertices;
  else
    return 0;
}

int
vtkFSSurfaceSource::FindVertexAtRAS ( float const iRAS[3], float* oDistance ) {

  float surf[3];
  this->ConvertRASToSurface( iRAS, surf );

  return this->FindVertexAtSurfaceRAS( surf, oDistance );
}

int
vtkFSSurfaceSource::FindVertexAtSurfaceRAS ( float const iSurfaceRAS[3],
					     float* oDistance ) {

  VERTEX v;
  v.x = iSurfaceRAS[0];
  v.y = iSurfaceRAS[1];
  v.z = iSurfaceRAS[2];
  float distance;
  int nClosestVertex =
    MHTfindClosestVertexNo( mHashTable, mMRIS, &v, &distance );

  if ( -1 == nClosestVertex ) {
    throw runtime_error( "No vertices found.");
  }

  if ( NULL != oDistance ) {
    *oDistance = distance;
  }

  return nClosestVertex;
}

void
vtkFSSurfaceSource::Execute () {

  if ( mMRIS == NULL) {
    return;
  }

  vtkPolyData *output = this->GetOutput();

  // Allocate all our arrays.
  int cVertices = mMRIS->nvertices;
  int cFaces = mMRIS->nfaces;
  vtkPoints* newPoints = vtkPoints::New();
  newPoints->Allocate( cVertices );
  vtkCellArray* newPolys = vtkCellArray::New();
  newPolys->Allocate( newPolys->EstimateSize(cFaces,VERTICES_PER_FACE) );
  vtkFloatArray* newNormals = vtkFloatArray::New();
  newNormals->Allocate( cVertices );
  newNormals->SetNumberOfComponents( 3 );
  newNormals->SetName( "Normals" );

  // Go through the surface and copy the vertex and normal for each
  // vertex.
  float point[3], normal[3];
  for ( int vno = 0; vno < cVertices; vno++ ) {

    point[0] = mMRIS->vertices[vno].x;
    point[1] = mMRIS->vertices[vno].y;
    point[2] = mMRIS->vertices[vno].z;
    newPoints->InsertNextPoint( point );

    normal[0] = mMRIS->vertices[vno].nx;
    normal[1] = mMRIS->vertices[vno].ny;
    normal[2] = mMRIS->vertices[vno].nz;
    newNormals->InsertNextTuple( normal );
  }

  // Go through and add the face indices.
  vtkIdType face[VERTICES_PER_FACE];
  for ( int fno = 0; fno < cFaces; fno++ ) {

    face[0] = mMRIS->faces[fno].v[0];
    face[1] = mMRIS->faces[fno].v[1];
    face[2] = mMRIS->faces[fno].v[2];
    newPolys->InsertNextCell( 3, face );
  }

  output->SetPoints( newPoints );
  newPoints->Delete();

  output->GetPointData()->SetNormals( newNormals );
  newNormals->Delete();

  newPolys->Squeeze(); // since we've estimated size; reclaim some space
  output->SetPolys( newPolys );
  newPolys->Delete();

}

