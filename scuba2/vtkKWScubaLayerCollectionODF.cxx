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
