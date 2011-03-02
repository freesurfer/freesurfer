/**
 * @file  vtkFSSurfaceSource.cxx
 * @brief import a freesurfer surface file into vtk
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:56 $
 *    $Revision: 1.14 $
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

#include <limits>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include "vtkFSSurfaceSource.h"
#include "vtkObjectFactory.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"

using namespace std;

vtkStandardNewMacro( vtkFSSurfaceSource );
vtkCxxRevisionMacro( vtkFSSurfaceSource, "$Revision: 1.14 $" );

vtkFSSurfaceSource::vtkFSSurfaceSource() :
    mMRIS( NULL ),
    mbBoundsCacheDirty( true ),
    mHashTable( NULL ) {

   this->vtkSource::SetNthOutput(0, vtkPolyData::New());
   
   // Releasing data for pipeline parallism.
   // Filters will know it is empty.
   this->Outputs[0]->ReleaseData();
   this->Outputs[0]->Delete();
}

vtkFSSurfaceSource::~vtkFSSurfaceSource () {

  if( NULL != mHashTable ) MHTfree( &mHashTable );
  if( NULL != mMRIS ) MRISfree( &mMRIS );
}

void
vtkFSSurfaceSource::MRISRead( char const* ifn ) {

  char* fn = strdup( ifn );
  MRIS* mris = ::MRISread( fn );
  free( fn );
  if ( mris == NULL ) {
    throw runtime_error( "MRISread failed" );
  }

  // Out with the old and in with the new.
  if( NULL != mMRIS ) {
    MRISfree( &mMRIS );
  }    
  mMRIS = mris;

  // Get some info from the MRIS. This can either come from the volume
  // geometry data embedded in the surface; this is done for newer
  // surfaces. Or it can come from the source information in the
  // transform. We use it to get the RAS center offset for the
  // surface->RAS transform.
  if ( mMRIS->vg.valid ) {

    mSurfaceToRASMatrix[0] = 1;
    mSurfaceToRASMatrix[1] = 0;
    mSurfaceToRASMatrix[2] = 0;
    mSurfaceToRASMatrix[3] = mMRIS->vg.c_r;
    mSurfaceToRASMatrix[4] = 0;
    mSurfaceToRASMatrix[5] = 1;
    mSurfaceToRASMatrix[6] = 0;
    mSurfaceToRASMatrix[7] = mMRIS->vg.c_a;
    mSurfaceToRASMatrix[8] = 0;
    mSurfaceToRASMatrix[9] = 0;
    mSurfaceToRASMatrix[10] = 1;
    mSurfaceToRASMatrix[11] = mMRIS->vg.c_s;
    mSurfaceToRASMatrix[12] = 0;
    mSurfaceToRASMatrix[13] = 0;
    mSurfaceToRASMatrix[14] = 0;
    mSurfaceToRASMatrix[15] = 1;
    
  } else if ( mMRIS->lta ) {

    mSurfaceToRASMatrix[0] = 1;
    mSurfaceToRASMatrix[1] = 0;
    mSurfaceToRASMatrix[2] = 0;
    mSurfaceToRASMatrix[3] = -mMRIS->lta->xforms[0].src.c_r;
    mSurfaceToRASMatrix[4] = 0;
    mSurfaceToRASMatrix[5] = 1;
    mSurfaceToRASMatrix[6] = 0;
    mSurfaceToRASMatrix[7] = -mMRIS->lta->xforms[0].src.c_a;
    mSurfaceToRASMatrix[8] = 0;
    mSurfaceToRASMatrix[9] = 0;
    mSurfaceToRASMatrix[10] = 1;
    mSurfaceToRASMatrix[11] = -mMRIS->lta->xforms[0].src.c_s;
    mSurfaceToRASMatrix[12] = 0;
    mSurfaceToRASMatrix[13] = 0;
    mSurfaceToRASMatrix[14] = 0;
    mSurfaceToRASMatrix[15] = 1;

  }

  // Make our transform object and set the matrix.
  mSurfaceToRASTransform = vtkSmartPointer<vtkTransform>::New();
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
vtkFSSurfaceSource::ConvertSurfaceToRAS ( double iX, double iY, double iZ,
					  double& oX, double& oY, double& oZ ) const {

  double surface[3];
  double ras[3];

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
vtkFSSurfaceSource::ConvertRASToSurface ( double iX, double iY, double iZ,
					  double& oX, double& oY, double& oZ ) const {

  double ras[3];
  double surface[3];
  
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
vtkFSSurfaceSource::ConvertSurfaceToRAS ( double const iSurf[3], double oRAS[3] ) const {

  mSurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void
vtkFSSurfaceSource::ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const {

  mSurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void
vtkFSSurfaceSource::ConvertRASToSurface ( double const iRAS[3], double oSurf[3] ) const {

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
vtkFSSurfaceSource::FindVertexAtRAS ( double const iRAS[3], double* oDistance ) {

  double surf[3];
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

int
vtkFSSurfaceSource::FindVertexAtSurfaceRAS ( double const iSurfaceRAS[3],
					     double* oDistance ) {

  VERTEX v;
  v.x = static_cast<float>(iSurfaceRAS[0]);
  v.y = static_cast<float>(iSurfaceRAS[1]);
  v.z = static_cast<float>(iSurfaceRAS[2]);
  float distance;
  int nClosestVertex =
    MHTfindClosestVertexNo( mHashTable, mMRIS, &v, &distance );

  if ( -1 == nClosestVertex ) {
    throw runtime_error( "No vertices found.");
  }

  if ( NULL != oDistance ) {
    *oDistance = static_cast<double>(distance);
  }

  return nClosestVertex;
}

void
vtkFSSurfaceSource::GetRASAtVertex ( int inVertex, float ioRAS[3] ) {

  float surfaceRAS[3];
  this->GetSurfaceRASAtVertex( inVertex, surfaceRAS );

  this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
}

void
vtkFSSurfaceSource::GetRASAtVertex ( int inVertex, double ioRAS[3] ) {

  double surfaceRAS[3];
  this->GetSurfaceRASAtVertex( inVertex, surfaceRAS );

  this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
}

void
vtkFSSurfaceSource::GetSurfaceRASAtVertex ( int inVertex, float ioRAS[3] ) {

  if( mMRIS == NULL )
    throw runtime_error( "GetRASAtVertex: mMRIS was NULL" );

  if( inVertex < 0 || inVertex >= mMRIS->nvertices )
    throw runtime_error( "GetRASAtVertex: inVertex was invalid" );

  ioRAS[0] = mMRIS->vertices[inVertex].x;
  ioRAS[1] = mMRIS->vertices[inVertex].y;
  ioRAS[2] = mMRIS->vertices[inVertex].z;
}

void
vtkFSSurfaceSource::GetSurfaceRASAtVertex ( int inVertex, double ioRAS[3] ) {

  if( mMRIS == NULL )
    throw runtime_error( "GetRASAtVertex: mMRIS was NULL" );

  if( inVertex < 0 || inVertex >= mMRIS->nvertices )
    throw runtime_error( "GetRASAtVertex: inVertex was invalid" );

  ioRAS[0] = static_cast<double>(mMRIS->vertices[inVertex].x);
  ioRAS[1] = static_cast<double>(mMRIS->vertices[inVertex].y);
  ioRAS[2] = static_cast<double>(mMRIS->vertices[inVertex].z);
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

  vtkSmartPointer<vtkPoints> newPoints = 
    vtkSmartPointer<vtkPoints>::New();
  newPoints->Allocate( cVertices );

  vtkSmartPointer<vtkCellArray> newPolys = 
    vtkSmartPointer<vtkCellArray>::New();
  newPolys->Allocate( newPolys->EstimateSize(cFaces,VERTICES_PER_FACE) );

  vtkSmartPointer<vtkFloatArray> newNormals =
    vtkSmartPointer<vtkFloatArray>::New();
  newNormals->Allocate( cVertices );
  newNormals->SetNumberOfComponents( 3 );
  newNormals->SetName( "Normals" );

  // Go through the surface and copy the vertex and normal for each
  // vertex. We have to transform them from surface RAS into normal
  // RAS.
  float point[3], normal[3], surfaceRAS[3];
  for ( int vno = 0; vno < cVertices; vno++ ) {

    surfaceRAS[0] = mMRIS->vertices[vno].x;
    surfaceRAS[1] = mMRIS->vertices[vno].y;
    surfaceRAS[2] = mMRIS->vertices[vno].z;
    this->ConvertSurfaceToRAS( surfaceRAS, point );
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
  output->GetPointData()->SetNormals( newNormals );
  newPolys->Squeeze(); // since we've estimated size; reclaim some space
  output->SetPolys( newPolys );
}

void
vtkFSSurfaceSource::FindPath ( int inStartVertex, int inEndVertex, 
			       vector<int>& iolPath ) {

  if( NULL == mMRIS )
    throw runtime_error( "Must have read in a surface before attempting to find a path" );

  const int cVerticies = mMRIS->nvertices;

  if( inStartVertex < 0 || inStartVertex >= cVerticies )
    throw runtime_error( "Start VNO is invalid" );
  if( inEndVertex < 0 || inEndVertex >= cVerticies )
    throw runtime_error( "End VNO is invalid" );

  // We use these arrays to hold the state of the search. The distance
  // is the distance from the start to that vertex to that vertex in
  // the current shortest path we've gound. The predecessor for a
  // vertex is the vertex number right before that vertex in the
  // shortest path. The check vector is a list of vertices to check
  // (our 'cloud').
  vector<float> aDistance( cVerticies, numeric_limits<float>::max() );
  vector<int> aPredecessor( cVerticies, -1 );
  vector<int> aCheck;
 
  // Start at the start vertex.
  aDistance[inStartVertex] = 0;
  aPredecessor[inStartVertex] = inStartVertex;
  aCheck.push_back( inStartVertex );

  // While we're not done and have things to check...
  bool bDone = false;
  vector<int>::iterator tCheck;
  float closestDistance = numeric_limits<float>::max();
  vector<int>::iterator tClosestVertex;
  int nClosestVertex;
  while( !bDone && !aCheck.empty() ) {
    
    // Find the closest vertex that needs checking.
    closestDistance = numeric_limits<float>::max();
    nClosestVertex = -1;
    for ( tCheck = aCheck.begin(); tCheck != aCheck.end(); ++tCheck ) {
      if( aDistance[*tCheck] < closestDistance ) {
	closestDistance = aDistance[*tCheck];
	tClosestVertex = tCheck;
      }
    }
    
    // Take out the closest vertex.
    nClosestVertex = *tClosestVertex;
    aCheck.erase( tClosestVertex );
    
    // If this is it, we're done!
    if( inEndVertex == nClosestVertex ) {

      bDone = true;

    } else {
      
      // Otherwise, look at all our neighbors. We'll call this vertex v.
      VERTEX* v = &(mMRIS->vertices[nClosestVertex]);
      for( int nNeighbor = 0; nNeighbor < v->vnum; nNeighbor++ ) {
	
	// Get a neighbor. We'll call it u.
	int nNeighborVertex = v->v[nNeighbor];
	VERTEX* u = &(mMRIS->vertices[nNeighborVertex]);
	
	// Calc the vector from u to v.
	float vuX = u->x - v->x;
	float vuY = u->y - v->y;
	float vuZ = u->z - v->z;

	// Calc the distance here.
	float distance = sqrt( vuX*vuX + vuY*vuY + vuZ*vuZ );

	// If this is a new shortest path to this vertex, update the
	// predecessor and the distance here, and add it to the list
	// of vertices to check next.
	if( distance + aDistance[nClosestVertex] < 
	    aDistance[nNeighborVertex] ) {
	  aPredecessor[nNeighborVertex] = nClosestVertex;
	  aDistance[nNeighborVertex] = distance + aDistance[nClosestVertex];
	  aCheck.push_back( nNeighborVertex );
	}
      }
    }
  }

  // Add the predecessors from the dest to the src in the output path.
  iolPath.clear();
  int nPathVertex = inEndVertex;
  iolPath.push_back( inEndVertex );
  while( aPredecessor[nPathVertex] != inStartVertex ) {
    
    iolPath.push_back( aPredecessor[nPathVertex] );
    nPathVertex = aPredecessor[nPathVertex];
  }
  iolPath.push_back( inStartVertex );

}

void
vtkFSSurfaceSource::SmoothValuesOnSurface ( vtkFloatArray& iValues,
					    int icSteps ) {

  assert( mMRIS );

  if( iValues.GetNumberOfTuples() != mMRIS->nvertices ) 
    throw runtime_error( "Number of tuples in values must be equal to the number of vertices in the surface." );

  if( icSteps <= 0 )
    throw runtime_error( "Number of steps must be > 0." );

  float* aTmpValues = new float[iValues.GetNumberOfTuples()];

  // This could probably done in a cleaner way with a VTK object, but
  // I want to use this code because it's exactly the same as the old
  // tksurfer method.
  for( int nStep = 0; nStep < icSteps; nStep++ ) {

    // Copy the current values into our temp array.
    for( int nValue = 0; nValue < iValues.GetNumberOfTuples(); nValue++ )
      aTmpValues[nValue] = iValues.GetTuple1( nValue );

    // For each vertex...
    for( int nVertex = 0; nVertex < mMRIS->nvertices; nVertex++ ) {

      // Sum up the values of this vertex and its neighbors.
      VERTEX* v = &mMRIS->vertices[nVertex];
      float sum = aTmpValues[nVertex];

      for( int nNeighbor = 0; nNeighbor < v->vnum; nNeighbor++ )
	sum += aTmpValues[v->v[nNeighbor]];

      // Set the average in the values array.
      iValues.SetTuple1( nVertex, sum / (float)(1 + v->vnum) );
    }
  }

  delete [] aTmpValues;
}
