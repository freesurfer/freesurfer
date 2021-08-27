/**
 * @brief Inflates (or deflates) a surface by moving verts along the normals
 *
 * This VTK filter class takes vtkPolyData as in put and outputs
 * vtkPolyData. It translates all vertices along the normal by
 * InflateFactor. A >0 value will inflate the poly data, and a <0
 * value will deflate it. Originally taken from vtkCleanPolyData.cxx.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "vtkInflatePolyData.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangle.h"

vtkStandardNewMacro(vtkInflatePolyData);

vtkInflatePolyData::vtkInflatePolyData() :
  InflateFactor( 0.1 ) {
}


vtkInflatePolyData::~vtkInflatePolyData() {
}


int
vtkInflatePolyData::RequestInformation ( vtkInformation *vtkNotUsed(iRequest),
			      vtkInformationVector **vtkNotUsed(ioaInputInfo),
				       vtkInformationVector *ioaOutputInfo ) {

  // Get the first info from the output. (Should only have 1.)
  vtkInformation *outInfo = ioaOutputInfo->GetInformationObject( 0 );

  // Just 1 piece.
  outInfo->
    Set( vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1 );
  
  return 1;
}


int
vtkInflatePolyData::RequestUpdateExtent ( vtkInformation *vtkNotUsed(iRequest),
					  vtkInformationVector **ioaInputInfo,
					vtkInformationVector *ioaOutputInfo ) {

  // Get the first info from input and output. (Should only have 1.)
  vtkInformation *inInfo = ioaInputInfo[0]->GetInformationObject( 0 );
  vtkInformation *outInfo = ioaOutputInfo->GetInformationObject( 0 );
  
  if( outInfo->
      Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) == 0 ) {

    inInfo->
      Set( vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), 0 );
    inInfo->
      Set( vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1 );

  } else {

    inInfo->
      Set( vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), -1 );
    inInfo->
      Set( vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 0 );
  }
  
  return 1;
}


int
vtkInflatePolyData::RequestData ( vtkInformation *vtkNotUsed(iRequest),
				  vtkInformationVector **ioaInputInfo,
				  vtkInformationVector *ioaOutputInfo ) {

  // Get the first info from input and output. (Should only have 1.)
  vtkInformation *inInfo = ioaInputInfo[0]->GetInformationObject( 0 );
  vtkInformation *outInfo = ioaOutputInfo->GetInformationObject( 0 );

  // Get the poly data from each.
  vtkPolyData* inputPolyData = 
    vtkPolyData::SafeDownCast( inInfo->Get(vtkDataObject::DATA_OBJECT()) );
  vtkPolyData* outputPolyData = 
    vtkPolyData::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()) );

  // Pass through most of our data.
  outputPolyData->CopyStructure( inputPolyData );
  outputPolyData->GetCellData()->PassData( inputPolyData->GetCellData() );
  //outputPolyData->GetPointData()->PassData( inputPolyData->GetPointData() );

  // Get our input points.
  vtkPoints* points = inputPolyData->GetPoints();
  if( (NULL == points) || (points->GetNumberOfPoints() < 1) ) {
    vtkDebugMacro(<<"No input points");
    return 1;
  }
  
  // Get our input polygons.
  vtkCellArray* polys = inputPolyData->GetPolys();
  if( (NULL == polys) || (polys->GetNumberOfCells() < 1) ) {
    vtkDebugMacro(<<"No input polys");
    return 1;
  }

  // Allocate our output points. This is the only thing we change in
  // this filter, so we'll make a new array of points, populate them,
  // and set them in the output. All the other output is copied from
  // the input.
  vtkPoints* inflatedPoints = vtkPoints::New();
  inflatedPoints->SetNumberOfPoints( points->GetNumberOfPoints() );

  // Keep track of the normals for each polygon.
  std::map<int,std::map<int,double> > aNormals;

  // For each polygon...
  vtkIdType cellID = 0;
  vtkIdType cPoints = 0;
  vtkIdType* pPoints = NULL;
  for( polys->InitTraversal(); 
       polys->GetNextCell( cPoints, pPoints ); cellID++ ) {
    
    // Compute the normal for this polygon.
    double normal[3] = { 0, 0, 0 };
    vtkPolygon::ComputeNormal( points, cPoints, pPoints, normal );
    vtkMath::Normalize( normal );
    
    // Store the normal.
    aNormals[cellID][0] = normal[0];
    aNormals[cellID][1] = normal[1];
    aNormals[cellID][2] = normal[2];
    
  }

  // For each point...
  vtkIdList* lCells = vtkIdList::New();
  lCells->Allocate( 3 );
  for( int pointID = 0; pointID < points->GetNumberOfPoints(); pointID++ ) {

    // Find the polygons that use this point.
    inputPolyData->GetPointCells( pointID, lCells );

    // Get an average normal for those polygons. For each cell, get
    // the normal we calculated earlier and add it to the sum.
    double sumNormal[3] = { 0, 0, 0 };
    for( int nCell = 0; nCell < lCells->GetNumberOfIds(); nCell++ ) {

      int cellID = lCells->GetId( nCell );
      sumNormal[0] += aNormals[cellID][0];
      sumNormal[1] += aNormals[cellID][1];
      sumNormal[2] += aNormals[cellID][2];
    }

    // Get the average normal.
    double averageNormal[3];
    averageNormal[0] = sumNormal[0] / lCells->GetNumberOfIds();
    averageNormal[1] = sumNormal[1] / lCells->GetNumberOfIds();
    averageNormal[2] = sumNormal[2] / lCells->GetNumberOfIds();
    vtkMath::Normalize( averageNormal );

    // Get this point.
    double point[3];
    points->GetPoint( pointID, point );

    // Move the point along that average normal by the inflate factor.
    point[0] += this->InflateFactor * averageNormal[0];
    point[1] += this->InflateFactor * averageNormal[1];
    point[2] += this->InflateFactor * averageNormal[2];

    // Set the translated point in the output points.
    inflatedPoints->SetPoint( pointID, point );

  }
  lCells->Delete();

  // Set the inflated points in the output.
  outputPolyData->SetPoints( inflatedPoints );
  inflatedPoints->Squeeze();
  inflatedPoints->Delete();

  // Pass along all the other data.
  outputPolyData->SetVerts( inputPolyData->GetVerts() );
  outputPolyData->SetLines( inputPolyData->GetLines() );
  outputPolyData->SetPolys( inputPolyData->GetPolys() );
  outputPolyData->SetStrips( inputPolyData->GetStrips() );

  this->Modified();

  return 1;
}

//--------------------------------------------------------------------------
void vtkInflatePolyData::PrintSelf(ostream& os, vtkIndent indent) 
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Inflation Factor: "
     << this->InflateFactor << endl;
}
