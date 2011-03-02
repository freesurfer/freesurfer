/**
 * @file  vtkKWScubaLayerCollectionDTI.cxx
 * @brief Implementation for DTI viewers.
 *
 */
/*
 * Original Author: Kevin Teich
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

#include "vtkKWScubaLayerCollectionDTI.h"

#include "vtkFSVolumeSource.h"
#include "vtkImageAppendComponents.h"
#include "vtkImageReslice.h"
#include "vtkKWCheckButton.h"
#include "vtkKWScubaLayer2DDTI.h"
#include "vtkKWScubaLayer3DDTI.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButtonSetWithLabel.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollectionDTI );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionDTI, "$Revision: 1.3 $" );

vtkKWScubaLayerCollectionDTI::vtkKWScubaLayerCollectionDTI () :
  mFAVolumeSource( NULL ),
  mEV1VolumeSource( NULL ),
  mEV2VolumeSource( NULL ),
  mEV3VolumeSource( NULL ),
  mEValuesVolumeSource( NULL ),
  mMergedSources( NULL ),
  mChkBtnRenderEdges( NULL ),
  mRadBtnScalingUniform( NULL ),
  mRadBtnScalingFA( NULL ),
  mRadBtnDetailLeast( NULL ),
  mRadBtnDetailLess( NULL ),
  mRadBtnDetailNormal( NULL ),
  mTensorInterpolation( VTK_RESLICE_LINEAR ),
  mbRenderEdges( false ),
  mScaling( ScubaCollectionPropertiesDTI::Uniform ),
  mDetail( ScubaCollectionPropertiesDTI::Less ),
  mfnFAVolume("") {
}

vtkKWScubaLayerCollectionDTI::~vtkKWScubaLayerCollectionDTI () {

}

void
vtkKWScubaLayerCollectionDTI::SetFAVolumeFileName ( const char* ifnVolume ) {

  // Set our file name and load the volume.
  mfnFAVolume = ifnVolume;
  this->LoadVolumesFromFileName();

  // Find a good label based on the filename and set it in the
  // layer.
  string fnVolume = ifnVolume;
  string::size_type lastSlash = fnVolume.rfind( "/" );
  this->SetLabel( fnVolume.substr( lastSlash+1, string::npos ).c_str() );

}

void
vtkKWScubaLayerCollectionDTI::AddControls ( vtkKWWidget* iPanel ) {

  // Tensor interpolation radio buttons 
  vtkKWRadioButtonSet* radBtnSetInterpolation = vtkKWRadioButtonSet::New();
  radBtnSetInterpolation->SetParent( iPanel );
  radBtnSetInterpolation->Create();
  radBtnSetInterpolation->PackHorizontallyOn();
  
  char sTclCmd[1024];

  vtkKWRadioButton* radBtnLinear = radBtnSetInterpolation->AddWidget( 0 );
  radBtnLinear->SetText( "Linear" );
  sprintf( sTclCmd, "SetTensorInterpolationType %d", VTK_RESLICE_LINEAR );
  radBtnLinear->SetCommand( this, sTclCmd );
  if ( mTensorInterpolation == VTK_RESLICE_LINEAR )
    radBtnLinear->SelectedStateOn();

  vtkKWRadioButton* radBtnCubic =
    radBtnSetInterpolation->AddWidget( 1 );
  radBtnCubic->SetText( "Cubic" );
  sprintf( sTclCmd, "SetTensorInterpolationType %d", VTK_RESLICE_CUBIC );
  radBtnCubic->SetCommand( this, sTclCmd );
  if ( mTensorInterpolation == VTK_RESLICE_CUBIC )
    radBtnCubic->SelectedStateOn();

  vtkKWRadioButton* radBtnNearestNeighbor =
    radBtnSetInterpolation->AddWidget( 2 );
  radBtnNearestNeighbor->SetText( "Nearest Neighbor" );
  sprintf( sTclCmd, "SetTensorInterpolationType %d", VTK_RESLICE_NEAREST );
  radBtnNearestNeighbor->SetCommand( this, sTclCmd );
  if ( mTensorInterpolation == VTK_RESLICE_NEAREST )
    radBtnNearestNeighbor->SelectedStateOn();
    
  radBtnLinear->Delete();
  radBtnCubic->Delete();
  radBtnNearestNeighbor->Delete();

  // Draw Edges check button --------------------------------------------
  mChkBtnRenderEdges = vtkKWCheckButton::New();
  mChkBtnRenderEdges->SetParent( iPanel );
  mChkBtnRenderEdges->Create();
  mChkBtnRenderEdges->SetAnchorToWest();
  mChkBtnRenderEdges->SetText( "Draw Edges" );
  mChkBtnRenderEdges->SetCommand( this, "SetRenderEdges" );
  if ( mbRenderEdges )
    mChkBtnRenderEdges->SelectedStateOn();
  // --------------------------------------------------------------------

  // Scaling radio buttons. ---------------------------------------------
  vtkKWRadioButtonSetWithLabel* radBtnSetTensorScaling =
    vtkKWRadioButtonSetWithLabel::New();
  radBtnSetTensorScaling->SetParent( iPanel );
  radBtnSetTensorScaling->Create();
  radBtnSetTensorScaling->GetWidget()->PackHorizontallyOn();
  radBtnSetTensorScaling->SetLabelText( "Scaling: " );

  mRadBtnScalingUniform = radBtnSetTensorScaling->GetWidget()->
    AddWidget( (int)ScubaCollectionPropertiesDTI::Uniform );
  mRadBtnScalingUniform->SetText( "Uniform" );
  sprintf( sTclCmd, "SetTensorScalingFromInt %d", 
	   (int)ScubaCollectionPropertiesDTI::Uniform );
  mRadBtnScalingUniform->SetCommand( this, sTclCmd );
  if ( mScaling == (int)ScubaCollectionPropertiesDTI::Uniform )
    mRadBtnScalingUniform->SelectedStateOn();

  mRadBtnScalingFA = radBtnSetTensorScaling->GetWidget()->
    AddWidget( (int)ScubaCollectionPropertiesDTI::FA );
  mRadBtnScalingFA->SetText( "FA" );
  sprintf( sTclCmd, "SetTensorScalingFromInt %d", 
	   (int)ScubaCollectionPropertiesDTI::FA );
  mRadBtnScalingFA->SetCommand( this, sTclCmd );
  if ( mScaling == (int)ScubaCollectionPropertiesDTI::FA )
    mRadBtnScalingFA->SelectedStateOn();
  // --------------------------------------------------------------------

  // Detail radio buttons. ---------------------------------------------
  vtkKWRadioButtonSetWithLabel* radBtnSetTensorDetail =
    vtkKWRadioButtonSetWithLabel::New();
  radBtnSetTensorDetail->SetParent( iPanel );
  radBtnSetTensorDetail->Create();
  radBtnSetTensorDetail->GetWidget()->PackHorizontallyOn();
  radBtnSetTensorDetail->SetLabelText( "Detail: " );

  mRadBtnDetailLeast = radBtnSetTensorDetail->GetWidget()->
    AddWidget( (int)ScubaCollectionPropertiesDTI::Least );
  mRadBtnDetailLeast->SetText( "Least" );
  sprintf( sTclCmd, "SetTensorDetailFromInt %d", 
	   (int)ScubaCollectionPropertiesDTI::Least );
  mRadBtnDetailLeast->SetCommand( this, sTclCmd );
  if ( mDetail == (int)ScubaCollectionPropertiesDTI::Least )
    mRadBtnDetailLeast->SelectedStateOn();

  mRadBtnDetailLess = radBtnSetTensorDetail->GetWidget()->
    AddWidget( (int)ScubaCollectionPropertiesDTI::Less );
  mRadBtnDetailLess->SetText( "Less" );
  sprintf( sTclCmd, "SetTensorDetailFromInt %d", 
	   (int)ScubaCollectionPropertiesDTI::Less );
  mRadBtnDetailLess->SetCommand( this, sTclCmd );
  if ( mDetail == (int)ScubaCollectionPropertiesDTI::Less )
    mRadBtnDetailLess->SelectedStateOn();

  mRadBtnDetailNormal = radBtnSetTensorDetail->GetWidget()->
    AddWidget( (int)ScubaCollectionPropertiesDTI::Normal );
  mRadBtnDetailNormal->SetText( "Normal" );
  sprintf( sTclCmd, "SetTensorDetailFromInt %d", 
	   (int)ScubaCollectionPropertiesDTI::Normal );
  mRadBtnDetailNormal->SetCommand( this, sTclCmd );
  if ( mDetail == (int)ScubaCollectionPropertiesDTI::Normal )
    mRadBtnDetailNormal->SelectedStateOn();

  // --------------------------------------------------------------------

  this->Script( "pack %s %s %s %s -side top -fill x",
                radBtnSetInterpolation->GetWidgetName(),
                mChkBtnRenderEdges->GetWidgetName(),
                radBtnSetTensorScaling->GetWidgetName(),
                radBtnSetTensorDetail->GetWidgetName() );

  radBtnSetTensorScaling->Delete();
  radBtnSetTensorDetail->Delete();
}

void
vtkKWScubaLayerCollectionDTI::RemoveControls () {

  if ( mChkBtnRenderEdges ) {
    mChkBtnRenderEdges->Delete();
    mChkBtnRenderEdges = NULL;
  }
  if ( mRadBtnScalingUniform ) {
    mRadBtnScalingUniform->Delete();
    mRadBtnScalingUniform = NULL;
  }
  if ( mRadBtnScalingFA ) {
    mRadBtnScalingFA->Delete();
    mRadBtnScalingFA = NULL;
  }
  if ( mRadBtnDetailLeast ) {
    mRadBtnDetailLeast->Delete();
    mRadBtnDetailLeast = NULL;
  }
  if ( mRadBtnDetailLess ) {
    mRadBtnDetailLess->Delete();
    mRadBtnDetailLess = NULL;
  }
  if ( mRadBtnDetailNormal ) {
    mRadBtnDetailNormal->Delete();
    mRadBtnDetailNormal = NULL;
  }
}

vtkImageAppendComponents* 
vtkKWScubaLayerCollectionDTI::GetMergedSource () const {
  return mMergedSources;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetFAVolumeSource() const {
  return mFAVolumeSource;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetEigenValueVolumeSource() const {
  return mEValuesVolumeSource;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetEigenVectorVolumeSource( const int iSource ) const {
  
  const int first = 1;
  const int second = 2;
  const int third = 3;
  
  vtkFSVolumeSource* source = NULL;
  
  if( first == iSource ) {
    source = this->GetEigenVector1VolumeSource();
  } else if( second == iSource ) {
    source = this->GetEigenVector2VolumeSource();    
  } else if( third == iSource ) {
    source = this->GetEigenVector3VolumeSource();    
  } 
  
  return source;
  
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetEigenVector1VolumeSource() const {
  return mEV1VolumeSource;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetEigenVector2VolumeSource() const {
  return mEV2VolumeSource;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionDTI::GetEigenVector3VolumeSource() const {
  return mEV3VolumeSource;
}

void
vtkKWScubaLayerCollectionDTI::SetRenderEdges ( int ibScaling ) {

  if( mbRenderEdges != ibScaling ) {
    mbRenderEdges = ibScaling;
    this->SendBroadcast( "RenderEdgesChanged", NULL );
  }
}

int 
vtkKWScubaLayerCollectionDTI::GetRenderEdges () const {
  return mbRenderEdges;
}

void
vtkKWScubaLayerCollectionDTI::SetTensorInterpolationType ( int iMode ) {
  if( mTensorInterpolation != iMode ) {
    mTensorInterpolation = iMode;
    this->SendBroadcast( "TensorInterpolationChanged", NULL );
  }
}

int 
vtkKWScubaLayerCollectionDTI::GetTensorInterpolationType () const {
  return mTensorInterpolation;
}

void 
vtkKWScubaLayerCollectionDTI::SetTensorScalingFromInt ( int iScaling ) {
  if( iScaling != ScubaCollectionPropertiesDTI::Uniform &&
      iScaling != ScubaCollectionPropertiesDTI::FA )
    throw runtime_error ( "SetTensorScalingTypeFromInt got invalid value" );
  
  this->SetTensorScaling( (ScubaCollectionPropertiesDTI::TensorScalingType)iScaling );
}

void 
vtkKWScubaLayerCollectionDTI::SetTensorScaling ( ScubaCollectionPropertiesDTI::TensorScalingType iScaling ) {
  if( mScaling != iScaling ) {
    mScaling = iScaling;
    this->SendBroadcast( "TensorScalingChanged" );
  }
}

ScubaCollectionPropertiesDTI::TensorScalingType 
vtkKWScubaLayerCollectionDTI::GetTensorScaling () const {
  return mScaling;
}


void 
vtkKWScubaLayerCollectionDTI::SetTensorDetailFromInt ( int iDetail ) {
  if( iDetail != ScubaCollectionPropertiesDTI::Least &&
      iDetail != ScubaCollectionPropertiesDTI::Less &&
      iDetail != ScubaCollectionPropertiesDTI::Normal )
    throw runtime_error( "SetTensorDetailTypeFromInt got invalid value" );
  
  this->SetTensorDetail( (ScubaCollectionPropertiesDTI::TensorDetailType) iDetail );
}

void 
vtkKWScubaLayerCollectionDTI::SetTensorDetail ( ScubaCollectionPropertiesDTI::TensorDetailType iDetail ) {
  if( mDetail != iDetail ) {
    mDetail = iDetail;
    this->SendBroadcast( "TensorDetailChanged" );
  }
}

ScubaCollectionPropertiesDTI::TensorDetailType 
vtkKWScubaLayerCollectionDTI::GetTensorDetail () const {
  return mDetail;
}



vtkKWScubaLayer*
vtkKWScubaLayerCollectionDTI::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {

  vtkKWScubaLayer* layer = NULL;
  if( vtkKWScubaView::TwoDee == iMode ) {

    vtkKWScubaLayer2DDTI* layer2DDTI = vtkKWScubaLayer2DDTI::New();
    layer2DDTI->SetDTIProperties( this );
    
    layer = (vtkKWScubaLayer*)layer2DDTI;

  } else if( vtkKWScubaView::ThreeDee == iMode ) {

    vtkKWScubaLayer3DDTI* layer3DDTI = vtkKWScubaLayer3DDTI::New();
    layer3DDTI->SetDTIProperties( this );
    
    layer = (vtkKWScubaLayer*)layer3DDTI;

  }

  return layer;
}

void
vtkKWScubaLayerCollectionDTI::LoadVolumesFromFileName () {

  // We have the FA file name in mfnFAVolume, but we'll need to figure
  // out the other ones.

  const char* sReplace =    "/fa";
  const char* sEV1Partial = "/eigvec1";
  const char* sEV2Partial = "/eigvec2";
  const char* sEV3Partial = "/eigvec3";
  const char* sEVPartial =  "/eigvals";

  string fnEV1 = mfnFAVolume;
  size_t nFound = fnEV1.find( sReplace );
  if( nFound != string::npos )
    fnEV1.replace( nFound, strlen( sReplace ), sEV1Partial );

  string fnEV2 = mfnFAVolume;
  nFound = fnEV2.find( sReplace );
  if( nFound != string::npos )
    fnEV2.replace( nFound, strlen( sReplace ), sEV2Partial );

  string fnEV3 = mfnFAVolume;
  nFound = fnEV3.find( sReplace );
  if( nFound != string::npos )
    fnEV3.replace( nFound, strlen( sReplace ), sEV3Partial );

  string fnEV = mfnFAVolume;
  nFound = fnEV.find( sReplace );
  if( nFound != string::npos )
    fnEV.replace( nFound, strlen( sReplace ), sEVPartial );

  mFAVolumeSource = vtkFSVolumeSource::New();
  mFAVolumeSource->MRIRead( mfnFAVolume.c_str() );
  mFAVolumeSource->Update();

  mEV1VolumeSource = vtkFSVolumeSource::New();
  mEV1VolumeSource->MRIRead( fnEV1.c_str() );
  mEV1VolumeSource->Update();
  
  mEV2VolumeSource = vtkFSVolumeSource::New();
  mEV2VolumeSource->MRIRead( fnEV2.c_str() );
  mEV2VolumeSource->Update();

  mEV3VolumeSource = vtkFSVolumeSource::New();
  mEV3VolumeSource->MRIRead( fnEV3.c_str() );
  mEV3VolumeSource->Update();

  mEValuesVolumeSource = vtkFSVolumeSource::New();
  mEValuesVolumeSource->MRIRead( fnEV.c_str() );
  mEValuesVolumeSource->Update();

  mMergedSources = vtkImageAppendComponents::New();
  mMergedSources->AddInputConnection( mFAVolumeSource->GetOutputPort() );  
  mMergedSources->AddInputConnection( mEV1VolumeSource->GetOutputPort() );
  mMergedSources->AddInputConnection( mEV2VolumeSource->GetOutputPort() );
  mMergedSources->AddInputConnection( mEV3VolumeSource->GetOutputPort() );
  mMergedSources->AddInputConnection( mEValuesVolumeSource->GetOutputPort() );

}
