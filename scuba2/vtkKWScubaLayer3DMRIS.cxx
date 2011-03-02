/**
 * @file  vtkKWScubaLayer3DMRIS.cxx
 * @brief Displays 3D surfaces
 *
 * Displays surfaces as a 3D model.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.7 $
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


#include <string>
#include <stdexcept>
#include "vtkKWScubaLayer3DMRIS.h"
#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesMRIS.h"
#include "ScubaInfoItem.h"
#include "vtkActor.h"
#include "vtkCutter.h"
#include "vtkDecimatePro.h"
#include "vtkFSSurfaceSource.h"
#include "vtkKWCheckButton.h"
#include "vtkKWLabel.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayer3DMRIS, "$Revision: 1.7 $" );

vtkKWScubaLayer3DMRIS::vtkKWScubaLayer3DMRIS () :
  mMRISProperties( NULL ),
  mbShowSurface( true ) {
}

vtkKWScubaLayer3DMRIS::~vtkKWScubaLayer3DMRIS () {
}

void
vtkKWScubaLayer3DMRIS::SetMRISProperties ( ScubaCollectionPropertiesMRIS const* iProperties ) {
  mMRISProperties = iProperties;
}


void
vtkKWScubaLayer3DMRIS::Create () {

}

void
vtkKWScubaLayer3DMRIS::LoadDataFromProperties () {

  // Bail if we don't have our source yet.
  if( NULL == mMRISProperties )
    throw runtime_error( "vtkKWScubaLayer3DMRIS::Create: No source" );
  
  //
  // Mappers for the 3D surface.
  //
  m3DNormalMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  m3DNormalMapper->SetInput( mMRISProperties->GetNormalModeOutput() );

  m3DFastMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  m3DFastMapper->SetInput( mMRISProperties->GetFastModeOutput() );

  if( mMRISProperties->GetScalarsValues() &&
      mMRISProperties->GetScalarsColors() ) {
    mMRISProperties->GetNormalModeOutput()->GetPointData()->
      SetScalars( reinterpret_cast<vtkDataArray*>(mMRISProperties->GetScalarsValues()) );
    m3DNormalMapper->SetLookupTable( mMRISProperties->GetScalarsColors() );
    m3DNormalMapper->UseLookupTableScalarRangeOn();
    mMRISProperties->GetFastModeOutput()->GetPointData()->
      SetScalars( reinterpret_cast<vtkDataArray*>(mMRISProperties->GetScalarsValues()) );
    m3DFastMapper->SetLookupTable( mMRISProperties->GetScalarsColors() );
    m3DFastMapper->UseLookupTableScalarRangeOn();
  }

  //
  // Actor for the 3D surface.
  //
  m3DActor = vtkSmartPointer<vtkActor>::New();
  m3DActor->SetMapper( m3DNormalMapper );
  m3DActor->SetProperty( m3DActor->MakeProperty() );
  m3DActor->GetProperty()->SetColor( 0.6, 0.5, 0.5 );

  // Add the prop if we're showing the surface.
  if( mbShowSurface )
    this->AddProp( m3DActor );

  // Update the fastClipper now so we do the decimation up front.
  m3DFastMapper->Update();


  float ras3DView[3];
  ras3DView[0] = mViewProperties->Get3DRASX();
  ras3DView[1] = mViewProperties->Get3DRASY();
  ras3DView[2] = mViewProperties->Get3DRASZ();

  // Now make three actors for the optional 2D plane intersection mode.
  for( int n = 0; n < 3; n++ ) {
    
    //
    // Cutting planes for the 2D intersections.
    //
    m2DSlicePlane[n] = vtkSmartPointer<vtkPlane>::New();
    m2DSlicePlane[n]->SetOrigin( (n==0) ? ras3DView[n] : 0,
				 (n==1) ? ras3DView[n] : 0,
				 (n==2) ? ras3DView[n] : 0 );
    m2DSlicePlane[n]->SetNormal( (n==0), (n==1), (n==2) );
    
    //
    // Cutters that takes the 3D surface and outputs 2D lines.
    //
    vtkSmartPointer<vtkCutter> clipper = 
      vtkSmartPointer<vtkCutter>::New();
    clipper->SetInputConnection( mMRISProperties->GetNormalModeOutputPort() );
    clipper->SetCutFunction( m2DSlicePlane[n] );
    
    vtkSmartPointer<vtkCutter> fastClipper =
      vtkSmartPointer<vtkCutter>::New();
    fastClipper->SetInputConnection( mMRISProperties->GetFastModeOutputPort() );
    fastClipper->SetCutFunction( m2DSlicePlane[n] );
    
    //
    // Mappers for the lines.
    //
    m2DNormalMapper[n] = vtkSmartPointer<vtkPolyDataMapper>::New();
    m2DNormalMapper[n]->SetInput( clipper->GetOutput() );
    
    m2DFastMapper[n] = vtkSmartPointer<vtkPolyDataMapper>::New();
    m2DFastMapper[n]->SetInput( fastClipper->GetOutput() );
    
    //
    // Actors in the scene, drawing the mapped lines.
    //
    m2DActor[n] = vtkSmartPointer<vtkActor>::New();
    m2DActor[n]->SetMapper( m2DNormalMapper[n] );
    m2DActor[n]->SetBackfaceProperty( m2DActor[n]->MakeProperty() );
    m2DActor[n]->GetBackfaceProperty()->BackfaceCullingOff();
    m2DActor[n]->SetProperty( m2DActor[n]->MakeProperty() );
    m2DActor[n]->GetProperty()->SetColor( 1, 0, 0 );
    m2DActor[n]->GetProperty()->SetEdgeColor( 1, 0, 0 );
    m2DActor[n]->GetProperty()->SetInterpolationToFlat();

    this->AddProp( m2DActor[n] );
  }
}

void
vtkKWScubaLayer3DMRIS::AddControls ( vtkKWWidget* iPanel ) {

  // Smooth check button ------------------------------------------------
  vtkSmartPointer<vtkKWCheckButton> chkBtnSurface = 
    vtkSmartPointer<vtkKWCheckButton>::New();
  chkBtnSurface->SetParent( iPanel );
  chkBtnSurface->Create();
  chkBtnSurface->SetAnchorToWest();
  chkBtnSurface->SetText( "Show Surface" );
  chkBtnSurface->SetCommand( this, "SetShowSurface" );
  if ( mbShowSurface )
    chkBtnSurface->SelectedStateOn();
  // --------------------------------------------------------------------

  this->Script( "pack %s -side top -fill x -anchor nw",
		chkBtnSurface->GetWidgetName() );
}

void
vtkKWScubaLayer3DMRIS::RemoveControls () {

}

void
vtkKWScubaLayer3DMRIS::DoListenToMessage ( string const isMessage,
					   void* const iData ) {
  
  if( isMessage == "OpacityChanged" ) {

    if ( m3DActor )
      if ( m3DActor->GetProperty() ) {
	m3DActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
	this->PipelineChanged();
      }
    
  } else if( isMessage == "FastModeChanged" ) {
    
//     if ( mProperties->GetFastMode() )
//       m3DActor->SetMapper( m3DFastMapper );
//     else
//       m3DNormalMapper->SetMapper( m3DNormalMapper );
    
//     this->PipelineChanged();
    
  } else if( isMessage == "Layer3DInfoChanged" ) {
    
    if ( m2DSlicePlane[0].GetPointer() &&
	 m2DSlicePlane[1].GetPointer() && 
	 m2DSlicePlane[2].GetPointer() ) {

      // For each plane, set the origin to the new RAS coord to draw
      // the new slice. Since the slice is drawn in the actual 3D
      // space of the cutting plane, we don't have to change anything
      // else. Sweet!
      float rasX = mViewProperties->Get3DRASX();
      m2DSlicePlane[0]->SetOrigin( rasX, 0, 0 );
      
      float rasY = mViewProperties->Get3DRASY();
      m2DSlicePlane[1]->SetOrigin( 0, rasY, 0 );
      
      float rasZ = mViewProperties->Get3DRASZ();
      m2DSlicePlane[2]->SetOrigin( 0, 0, rasZ );
      
      this->PipelineChanged();
    }
  }
}

void
vtkKWScubaLayer3DMRIS::GetRASBounds ( float ioBounds[6] ) const {

  if ( mMRISProperties )
    mMRISProperties->GetSource()->GetRASBounds( ioBounds );
  else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer3DMRIS::GetInfoItems ( float iRAS[3],
                                      list<ScubaInfoItem>& ilInfo ) const {

  ScubaInfoItem info;
  char sLabel[1024], sValue[1024];

  try {
    float distance;
    int nVertex =
      mMRISProperties->GetSource()->FindVertexAtSurfaceRAS( iRAS, &distance );

    sprintf( sLabel, "%s vno", mProperties->GetLabel() );
    sprintf( sValue, "%d (%.2f)", nVertex, distance );
    info.Clear();
    info.SetLabel( sLabel );
    info.SetValue( sValue );
    info.SetShortenHint( false );
  } catch (...) {

    // No vertex found, don't worry.
    sprintf( sLabel, "%s vno", mProperties->GetLabel() );
    strncpy( sValue, "NONE", sizeof(sValue) );
    info.Clear();
    info.SetLabel( sLabel );
    info.SetValue( sValue );
    info.SetShortenHint( false );
  }

  // Return it.
  ilInfo.push_back( info );
}

void
vtkKWScubaLayer3DMRIS::SetShowSurface ( int ibShow ) {

  if( ibShow != mbShowSurface ) {

    mbShowSurface = ibShow;

    if( mbShowSurface ) {
      this->AddProp( m3DActor );
    } else {
      this->RemoveProp( m3DActor );
    }

    this->PipelineChanged();
  }
}
