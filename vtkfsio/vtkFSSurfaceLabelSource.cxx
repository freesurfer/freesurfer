/**
 * @file  vtkFSSurfaceLabelSource.h
 * @brief Reads a label, maps it to a surface, and outputs PolyData
 *
 * A FreeSurfer label file consists of a list of vertices that may
 * also have associated vertex indices. This will read in a label, map
 * it to a surface, fill in any holes, and output PolyData, which will
 * appear to be a subset of the surface PolyData.
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

#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <assert.h>
#include "vtkFSSurfaceLabelSource.h"
#include "vtkFSSurfaceSource.h"
#include "vtkObjectFactory.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

using namespace std;

vtkStandardNewMacro( vtkFSSurfaceLabelSource );
vtkCxxRevisionMacro( vtkFSSurfaceLabelSource, "$Revision: 1.14 $" );

vtkFSSurfaceLabelSource::vtkFSSurfaceLabelSource() :
  LabelFileName( NULL ), Mris( NULL ), Label( NULL ) {

  this->vtkSource::SetNthOutput(0, vtkPolyData::New());

  // Releasing data for pipeline parallism.  Filters will know it is
  // empty.
  this->Outputs[0]->ReleaseData();
  this->Outputs[0]->Delete();
}

vtkFSSurfaceLabelSource::~vtkFSSurfaceLabelSource () {

}

vtkPolyData*
vtkFSSurfaceLabelSource::GetOutput () {

  if ( this->NumberOfOutputs < 1 ) {
    return NULL;
  }

  return (vtkPolyData *)(this->Outputs[0]);
}

vtkPolyData*
vtkFSSurfaceLabelSource::GetOutput ( int inOutput ) {

  return (vtkPolyData *)this->vtkSource::GetOutput( inOutput );
}

void
vtkFSSurfaceLabelSource::SetOutput ( vtkPolyData* iOutput ) {

  this->vtkSource::SetNthOutput( 0, iOutput );
}

void
vtkFSSurfaceLabelSource::InitializeEmptyLabel () {

  // Delete existing label if present, and create a new one with 0
  // vertices.
  if( NULL != Label ) {
    LabelFree( &Label );
  }

  Label = LabelAlloc( 0, (char*)"", (char*)"" );

  this->Modified();
}

void
vtkFSSurfaceLabelSource::AddVerticesToLabel ( int icVertices, 
					      int* iaVertices ) {
  assert( icVertices > 0 );
  assert( iaVertices );
  assert( Mris );
  assert( Label );

  // Reallocate the label with enough points.
  LabelRealloc( Label, Label->n_points + icVertices );

  // For each vertex to add...
  int nPoint = Label->n_points;
  for( int nVertex = 0; nVertex < icVertices; nVertex++ ) {
    
    // Check the vertex number.
    if( iaVertices[nVertex] < 0 ||
	iaVertices[nVertex] >= Mris->nvertices )
      throw runtime_error( "Invalid vertex number" );

    // Make a new point.
    Label->lv[nPoint].vno = iaVertices[nVertex];
    Label->lv[nPoint].x = Mris->vertices[iaVertices[nVertex]].x;
    Label->lv[nPoint].y = Mris->vertices[iaVertices[nVertex]].y;
    Label->lv[nPoint].z = Mris->vertices[iaVertices[nVertex]].z;
    Label->lv[nPoint].deleted = 0;

    Label->n_points++;
    nPoint++;
  }

  this->Modified();
}

void
vtkFSSurfaceLabelSource::RemoveVerticesFromLabel ( int icVertices, 
						   int* iaVertices ) {

  assert( Label );
  assert( iaVertices );

  // For each vertex to remove...
  for( int nVertex = 0; nVertex < icVertices; nVertex++ ) {
    
    // Check the vertex.
    if( iaVertices[nVertex] < 0 ||
	iaVertices[nVertex] >= Mris->nvertices )
      throw runtime_error( "Invalid vertex number" );

    // Search for this vertex in the label and if found, set its
    // deleted flag to true.
    for( int nPoint = 0; nPoint < Label->n_points; nPoint++ ) {
      
      if( Label->lv[nPoint].vno == iaVertices[nVertex] )
	Label->lv[nPoint].deleted = 1;

    }
  }

  this->Modified();
}

void
vtkFSSurfaceLabelSource::Execute () {
  
  if( NULL == Mris )
    vtkErrorMacro( << "vtkFSSurfaceLabelSource cannot exectue without a surface" );

  if( NULL == Label ) {
    this->ReadLabelFile();
    if( NULL == Label )
      vtkErrorMacro( << "vtkFSSurfaceLabelSource cannot exectue without a label" );
  }
  
  // We use marks to figure out where the label is on the surface.
  MRISclearMarks( Mris );

  LabelMarkUndeleted( Label, Mris );
  
  // Now we have a marked surface. What we want to do is create a
  // new poly data object without only marked faces in it. To do
  // that, we need a list of points in the faces, and make faces
  // that index those points. However we need to map the surface
  // based indices to a new range of indices that only have the
  // selected points. So we'll save a list of points and faces that
  // are selected using the old indices, and then create a mapping.
  
  // For each face in the surface, if all its vertices are marked,
  // put those points and the face in our list of ones that will go
  // into the new poly data.
  map<vtkIdType,map<int,double> > aPoints;
  map<vtkIdType,map<int,vtkIdType> > aPolys;
  int nPoly = 0;
  for( int nFace = 0; nFace < Mris->nfaces; nFace++ ) {
    if( Mris->vertices[Mris->faces[nFace].v[0]].marked &&
	Mris->vertices[Mris->faces[nFace].v[1]].marked &&
	Mris->vertices[Mris->faces[nFace].v[2]].marked ) {
      
      // We're going to make an entry for a poly with a 0-based face
      // index, but using the same vertex index as the original
      // surface.
      for( int nFaceVertex = 0; 
	   nFaceVertex < VERTICES_PER_FACE; nFaceVertex++ ) {
	
	int nVertex = Mris->faces[nFace].v[nFaceVertex];
	
	// Save the poly using the 0 based index.
	aPolys[nPoly][nFaceVertex] = nVertex;
	
	// Save the point with the original index. This can be
	// redundant, we don't care.
	aPoints[nVertex][0] = Mris->vertices[nVertex].x;
	aPoints[nVertex][1] = Mris->vertices[nVertex].y;
	aPoints[nVertex][2] = Mris->vertices[nVertex].z;
      }
      
      nPoly++;
    }
  }
  
  // Now lets make point and poly lists out of the maps we just
  // made.
  vtkPoints* labelPoints = vtkPoints::New();
  labelPoints->SetNumberOfPoints( aPoints.size() );
  
  vtkCellArray* labelPolys = vtkCellArray::New();
  labelPolys->Allocate( aPolys.size() );
  
  vtkIdType nNextNewID = 0;
  map<vtkIdType,vtkIdType> aOldIDToNewID;
  
  // For each point we saved, we need to map its surface based index
  // to a new 0 based index.
  map<vtkIdType,map<int,double> >::iterator tPoint;
  for( tPoint = aPoints.begin(); tPoint != aPoints.end(); ++tPoint ) {
    
    // Get the old ID.
    vtkIdType oldID = tPoint->first;
    
    // Build a point.
    double point[3];
    for( int n = 0; n < 3; n++ )
      point[n] = tPoint->second[n];
    
    // Map the old ID to a new ID and save it.
    aOldIDToNewID[oldID] = nNextNewID++;
    
    // Insert the point with the new ID.
    labelPoints->SetPoint( aOldIDToNewID[oldID], point );
  }
  
  // Now for each poly, add it to the polys array using the new 0
  // based point indices.
  map<vtkIdType,map<int,vtkIdType> >::iterator tPoly;
  for( tPoly = aPolys.begin(); tPoly != aPolys.end(); ++tPoly ) {
    
    // Make a poly using the new IDs.
    vtkIdType poly[3];
    for( int n = 0; n < 3; n++ )
      poly[n] = aOldIDToNewID[tPoly->second[n]];
    
    // Insert the poly.
    labelPolys->InsertNextCell( 3, poly );
  }

  // Set the points and polys in our output
  this->GetOutput()->SetPoints( labelPoints );
  labelPoints->Delete();
  
  this->GetOutput()->SetPolys( labelPolys );
  labelPolys->Delete();
  
}

void
vtkFSSurfaceLabelSource::ReadLabelFile () {

  if( NULL == LabelFileName )
    vtkErrorMacro( << "vtkFSSurfaceLabelSource cannot exectue without a label file name" );


  // Load the white vertex positions in the surface.
  int eMRIS = MRISreadWhiteCoordinates( Mris, (char*)"white" );
  if( eMRIS != 0 )
    throw runtime_error( "Couldn't read the white surface file, so unable to load a label." );
  
  // Load the label file.
  char* fnLabel = strdup( LabelFileName );
  Label = LabelRead( NULL, fnLabel );
  free( fnLabel );
  if( NULL == Label )
    throw runtime_error( "Couldn't read label" );
    
  // Map it to the surface using the white coordinates.
  LabelToWhite( Label, Mris );
  
  // See if it's completely unassigned. If so, there are likely to
  // be lots of holes in the label after we sample it to the
  // vertices, so we'll fill it in afterwards.
  int unassigned;
  LabelIsCompletelyUnassigned( Label, &unassigned );
  
  // Assign the mris vertex numbers to unnumbered vertices in the
  // label.
  LabelFillUnassignedVertices( Mris, Label, WHITE_VERTICES );
  
  // If we were unassigned before, fill the holes now.
  if( unassigned ) {
    LABEL* filledLabel = LabelFillHoles( Label, Mris, WHITE_VERTICES );
    LabelFree( &Label );
    Label = filledLabel;
  }
    
}

void
vtkFSSurfaceLabelSource::WriteLabelFile () {

  assert( Label );
  if( NULL == LabelFileName )
    vtkErrorMacro( << "vtkFSSurfaceLabelSource cannot write without a label file name" );

  // Write the file.
  char* fn = strdup( LabelFileName );
  int rLabel = LabelWrite( Label, fn );
  free( fn );

  if( rLabel != 0 ) {
    stringstream ssError;
    ssError << "Couldn't write the label to " << LabelFileName;
    throw runtime_error( ssError.str().c_str() );
  }
}

void
vtkFSSurfaceLabelSource::GetLabeledPoints ( vtkPoints& ioPoints ) {

  assert( Label );

  // Reset the point list and set the type to point.
  ioPoints.Reset();
  ioPoints.SetDataTypeToFloat();

  // Go through our Label, and for each not-deleted point, add it to
  // our points.
  for( int nPoint = 0; nPoint < Label->n_points; nPoint++ ) {
    
    if( !Label->lv[nPoint].deleted ) {
      
      float x = Label->lv[nPoint].x;
      float y = Label->lv[nPoint].y;
      float z = Label->lv[nPoint].z;
    
      ioPoints.InsertNextPoint( x, y, z );
    }
  }
}

void
vtkFSSurfaceLabelSource::GetLabeledVertices ( vector<int>& iolVertices ) {

  assert( Label );

  // Reset the list.
  iolVertices.clear();

  // Go through our Label, and for each not-deleted point, add the
  // vertex number to the list.
  for( int nPoint = 0; nPoint < Label->n_points; nPoint++ )
    if( !Label->lv[nPoint].deleted )
      iolVertices.push_back( Label->lv[nPoint].vno );

}
