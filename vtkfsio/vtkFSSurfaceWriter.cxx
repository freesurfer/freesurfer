/**
 * @file  vtkFSSurfaceWriter.h
 * @brief Writes a vtkPolyData to an MRIS file.
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/01/31 22:20:47 $
 *    $Revision: 1.3 $
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

#include <cassert>

#include "vtkFSSurfaceWriter.h"

#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

extern "C" {
#include "mrisurf.h"
}

vtkCxxRevisionMacro(vtkFSSurfaceWriter, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkFSSurfaceWriter);

vtkFSSurfaceWriter::vtkFSSurfaceWriter() :
  mMRIS( NULL ),
  FileName( NULL ) {
}

vtkFSSurfaceWriter::~vtkFSSurfaceWriter() {

  if( NULL != mMRIS )
    MRISfree( &mMRIS );
}

MRIS*
vtkFSSurfaceWriter::GetMRIS () {

  return mMRIS;
}

void
vtkFSSurfaceWriter::WriteData () {

  vtkPolyData* input = this->GetInput();
  if( NULL == input )
    return;
  
  // Try to allocate a surface.
  MRIS* mris = 
    MRISalloc( input->GetNumberOfPoints(), input->GetNumberOfPolys() );
  if( NULL == mris ) {
    vtkErrorMacro("Could not allocate MRIS with "
		  << input->GetNumberOfPoints() << " points and "
		  << input->GetNumberOfPolys() << " faces." );
    return;
  }

  // Make sure we write the proper kind of surface.
  mris->type = MRIS_TRIANGULAR_SURFACE;

  // Copy in the vertices.
  vtkPoints* points = input->GetPoints();
  assert( points );
  int cPoints = input->GetNumberOfPoints();
  for( int nPoint = 0; nPoint < cPoints; nPoint++ ) {

    double* point = points->GetPoint( nPoint );
    mris->vertices[nPoint].x = point[0];
    mris->vertices[nPoint].y = point[1];
    mris->vertices[nPoint].z = point[2];
  }
  
  // Copy in the faces.
  vtkIdType nFace = 0;
  vtkIdType cPointIDs = 0;
  vtkIdType* pPointIDs = NULL;
  vtkCellArray* polys = input->GetPolys();
  assert( polys );
  for( polys->InitTraversal(); 
       polys->GetNextCell( cPointIDs, pPointIDs ); nFace++ ) {

    if( cPointIDs != 3 ) {
      vtkErrorMacro("Face with invalid number of points: face "
		    << nFace << " with " << cPointIDs << " points." );
      MRISfree( &mris );
      return;
    }

    for( int nPointID = 0; nPointID < 3; nPointID++ )
      mris->faces[nFace].v[nPointID] = pPointIDs[nPointID];
  }

  // Delete the last MRIS and save this one.
  if( NULL != mMRIS )
    MRISfree( &mMRIS );
  mMRIS = mris;

  // Write the data.
  char* fnMRIS = strdup( this->FileName );
  MRISwrite( mMRIS, fnMRIS );
  free( fnMRIS );
}

int
vtkFSSurfaceWriter::FillInputPortInformation ( int, vtkInformation* ioInfo ) {

  ioInfo->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData" );
  return 1;
}

vtkPolyData*
vtkFSSurfaceWriter::GetInput () {
 
  return vtkPolyData::SafeDownCast( this->Superclass::GetInput() );
}

vtkPolyData*
vtkFSSurfaceWriter::GetInput ( int iPort ) {

  return vtkPolyData::SafeDownCast( this->Superclass::GetInput( iPort ) );
}

void
vtkFSSurfaceWriter::PrintSelf ( ostream& iStream, vtkIndent iIndent ) {

  this->Superclass::PrintSelf( iStream, iIndent );
}
