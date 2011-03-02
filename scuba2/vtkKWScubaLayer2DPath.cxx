/**
 * @file  vtkKWScubaLayer2DPath.cxx
 * @brief Displays 2D path
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.3 $
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


#include "vtkKWScubaLayer2DPath.h"

#include <string>
#include <stdexcept>

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesPath.h"
#include "ScubaInfoItem.h"

#include "vtkFSVolumeSource.h"

#include "vtkActor.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkContourFilter.h"
#include "vtkProperty.h"
#include "vtkPolyDataNormals.h"
#include "vtkStripper.h"
#include "vtkImageReslice.h"
#include "vtkTransform.h"
#include "vtkPlane.h"
#include "vtkCutter.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer2DPath );
vtkCxxRevisionMacro( vtkKWScubaLayer2DPath, "$Revision: 1.3 $" );

vtkKWScubaLayer2DPath::vtkKWScubaLayer2DPath () :
  mPathProperties( NULL ),
  mTriangleStripper( NULL ),
  mNormals( NULL ) {
};

vtkKWScubaLayer2DPath::~vtkKWScubaLayer2DPath () {
  if( mTriangleStripper ){
    mTriangleStripper->Delete();
  }
  
}

void
vtkKWScubaLayer2DPath::SetPathProperties ( ScubaCollectionPropertiesPath const* iProperties ) {
  mPathProperties = iProperties;
}


void
vtkKWScubaLayer2DPath::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mPathProperties )
    throw runtime_error( "vtkKWScubaLayer2DPath::Create: No source" );

  // Get some info from the MRIS.
  mRASCenter[0] = mPathProperties->GetPathVolumeSource()->GetRASCenterX();
  mRASCenter[1] = mPathProperties->GetPathVolumeSource()->GetRASCenterY();
  mRASCenter[2] = mPathProperties->GetPathVolumeSource()->GetRASCenterZ();

  // create triangle strips for faster rendering
  mTriangleStripper = vtkStripper::New();
  mTriangleStripper->SetInput( mPathProperties->GetMesh() );

  mNormals = vtkPolyDataNormals::New();
  mNormals->SetInputConnection( mTriangleStripper->GetOutputPort() );
  
  //
  // Cutting plane.
  //
  const float rasZ = mViewProperties->Get2DRASZ();
  const int inPlane = mViewProperties->Get2DInPlane();
  mSlicePlane = vtkPlane::New();
  mSlicePlane->SetOrigin( rasZ, rasZ, rasZ );
  mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );

  //
  // Cutters that takes the 3D surface and outputs 2D lines.
  //
  vtkCutter* clipper = vtkCutter::New();
  clipper->SetInputConnection( mNormals->GetOutputPort() );
  clipper->SetCutFunction( mSlicePlane );

  //
  // Mappers for the lines.
  //
  mMapper = vtkPolyDataMapper::New();
  mMapper->SetInput( clipper->GetOutput() );
  mMapper->ScalarVisibilityOff();
  

  //
  // Actors in the scene, drawing the mapped lines.
  //
  mActor = vtkActor::New();
  mActor->SetMapper( mMapper );
  mActor->SetBackfaceProperty( mActor->MakeProperty() );
  mActor->GetBackfaceProperty()->BackfaceCullingOff();
  mActor->SetProperty( mActor->MakeProperty() );
  mActor->GetProperty()->SetInterpolationToFlat();

  // get the color from the collection
  this->UpdatePathColor();

  this->AddProp( mActor );
      
}

void
vtkKWScubaLayer2DPath::GetRASBounds ( float ioBounds[6] ) const {

  if ( mPathProperties ) {
    mPathProperties->GetPathVolumeSource()->GetRASBounds( ioBounds );
  } else {
    for ( int nBound = 0; nBound < 6; nBound++ ) {
      ioBounds[nBound] = 0;
    }
  }
  
}

void
vtkKWScubaLayer2DPath::Get2DRASZIncrementHint ( float ioHint[3]) const {

  ioHint[0] = 1.0;
  ioHint[1] = 1.0;
  ioHint[2] = 1.0;
}

void
vtkKWScubaLayer2DPath::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "OpacityChanged" ) {

    if ( mActor ) {
      if ( mActor->GetProperty() ) {
        mActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
        this->PipelineChanged();
      }
    }
    
  } else if( isMessage == "Layer2DInfoChanged" ) {

    if ( mSlicePlane ) {
      
      const float rasZ = mViewProperties->Get2DRASZ();
      const int inPlane = mViewProperties->Get2DInPlane();
    
      // Update the location of the slice plane. We also transform the
      // RAS by the world center here.
      mSlicePlane->SetOrigin( (inPlane==0) ? rasZ - mRASCenter[0] : 0,
            (inPlane==1) ? rasZ - mRASCenter[1] : 0,
            (inPlane==2) ? rasZ - mRASCenter[2] : 0);
      mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );
      
      // When we slice the polydata, the physical location of the lines
      // in 3D space is that of the RAS plane on which we sliced, so we
      // need to transform the resulting data onto the 0 plane. However,
      // for some weird reason, the surface will disappear under the
      // texture planes in some slices, so move it just a little in
      // front as well. Also adjust by the world center.
      mActor->SetPosition( (inPlane==0) ? -(rasZ - mRASCenter[0]) + 0.1 :
         mRASCenter[0],
         (inPlane==1) ? -(rasZ - mRASCenter[1]) + 0.1 :
         mRASCenter[1],
         (inPlane==2) ? -(rasZ - mRASCenter[2]) - 0.1 :
         mRASCenter[2] );
      
      this->PipelineChanged();
    }
  } else if( isMessage == "PathColorChanged" ) {
    
    this->UpdatePathColor();
    
  }
  
}

void 
vtkKWScubaLayer2DPath::UpdatePathColor () {
  if( NULL != mPathProperties ) {
    
    double r, g, b;
    mPathProperties->GetPathColor( r, g, b );

    if( mActor ) {
      mActor->GetProperty()->SetColor( r, g, b );
      mActor->GetProperty()->SetEdgeColor( r, g, b );
    }
  }
}


