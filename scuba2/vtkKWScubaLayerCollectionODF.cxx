/**
 * @file  vtkKWScubaLayerCollectionODF.cxx
 * @brief Implementation for ODF viewers.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *    $Revision$
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

#include "vtkKWScubaLayerCollectionODF.h"

#include "vtkFSVolumeSource.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollectionODF );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionODF, "$Revision$" );

vtkKWScubaLayerCollectionODF::vtkKWScubaLayerCollectionODF () :
  mODFVolumeSource( NULL ),
  mfnODFVolume("") {
}

vtkKWScubaLayerCollectionODF::~vtkKWScubaLayerCollectionODF () {
}

void
vtkKWScubaLayerCollectionODF::AddControls ( vtkKWWidget* iPanel ) {
}

void
vtkKWScubaLayerCollectionODF::RemoveControls () {
}

void 
vtkKWScubaLayerCollectionODF::SetODFVolumeFileName ( const char* ifnVolume ) {
  mfnODFVolume = ifnVolume;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionODF::GetODFVolumeSource() const {
  return mODFVolumeSource;
}

vtkKWScubaLayer*
vtkKWScubaLayerCollectionODF::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {

  vtkKWScubaLayer* layer = NULL;

  if( vtkKWScubaView::TwoDee == iMode ) {

//    vtkKWScubaLayer2DDTI* layer2DDTI = vtkKWScubaLayer2DDTI::New();
//    layer2DDTI->SetDTIProperties( this );
//    
//    layer = (vtkKWScubaLayer*)layer2DDTI;

  } else if( vtkKWScubaView::ThreeDee == iMode ) {

//    vtkKWScubaLayer3DDTI* layer3DDTI = vtkKWScubaLayer3DDTI::New();
//    layer3DDTI->SetDTIProperties( this );
//    
//    layer = (vtkKWScubaLayer*)layer3DDTI;

  }

  return layer;
}

void
vtkKWScubaLayerCollectionODF::LoadVolumesFromFileName () {
}
