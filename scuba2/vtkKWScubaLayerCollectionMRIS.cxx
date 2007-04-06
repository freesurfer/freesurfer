/**
 * @file  vtkKWScubaLayerCollectionMRIS.cxx
 * @brief Implementation for MRIS viewers.
 *
 * In 2D, the MRIS is viewed as the intersection of the surface with a
 * slice. In 3D, the MRIS is viewed as a fully 3D mesh surface.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
 *    $Revision: 1.1 $
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

#include "vtkKWScubaLayerCollectionMRIS.h"
#include "vtkKWScubaLayer2DMRIS.h"
#include "vtkKWScubaLayer3DMRIS.h"
#include "vtkObjectFactory.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkDecimatePro.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkCommand.h"
#include "vtkKWProgressDialog.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollectionMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionMRIS, "$Revision: 1.1 $" );

vtkKWScubaLayerCollectionMRIS::vtkKWScubaLayerCollectionMRIS () :
  mSource( NULL ),
  mDecimator( NULL ),
  mfnSurface("") {

}

vtkKWScubaLayerCollectionMRIS::~vtkKWScubaLayerCollectionMRIS () {

}

void
vtkKWScubaLayerCollectionMRIS::SetSurfaceFileName ( const char* ifnSurface ) {

  // Set our file name and load the surface.
  mfnSurface = ifnSurface;
  this->LoadSurfaceFromFileName();

  // Find a good label based on the filename and set it in the
  // layer.
  string fnSurface = ifnSurface;
  string::size_type lastSlash = fnSurface.rfind( "/" );
  this->SetLabel( fnSurface.substr( lastSlash+1, string::npos ).c_str() );

  //
  // Decimator.
  //
  vtkKWProgressDialog* d = vtkKWProgressDialog::New();
  d->SetApplication( this->GetApplication() );
  d->SetWindowTitle( "Decimating Surface" );
  mDecimator = vtkDecimatePro::New();
  mDecimator->AddObserver( vtkCommand::ProgressEvent, d );
  mDecimator->SetInputConnection( mSource->GetOutputPort() );
  mDecimator->SetTargetReduction( 0.9 );

  // Force the decimator to do its thing now.
  mDecimator->Update();
}

vtkFSSurfaceSource*
vtkKWScubaLayerCollectionMRIS::GetSource () const {
  return mSource;
}

vtkAlgorithmOutput*
vtkKWScubaLayerCollectionMRIS::GetNormalModeOutputPort () const {
  return mSource->GetOutputPort();
}

vtkAlgorithmOutput*
vtkKWScubaLayerCollectionMRIS::GetFastModeOutputPort () const {
  return mDecimator->GetOutputPort();
}

vtkPolyData*
vtkKWScubaLayerCollectionMRIS::GetNormalModeOutput () const {
  return mSource->GetOutput();
}

vtkPolyData*
vtkKWScubaLayerCollectionMRIS::GetFastModeOutput () const {
  return mDecimator->GetOutput();
}

vtkKWScubaLayer*
vtkKWScubaLayerCollectionMRIS::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {

  vtkKWScubaLayer* layer = NULL;
  if( vtkKWScubaView::TwoDee == iMode ) {
    
    vtkKWScubaLayer2DMRIS* layer2DMRIS = vtkKWScubaLayer2DMRIS::New();
    layer2DMRIS->SetMRISProperties( this );

    layer = (vtkKWScubaLayer*)layer2DMRIS;

  } else if( vtkKWScubaView::ThreeDee == iMode ) {

    vtkKWScubaLayer3DMRIS* layer3DMRIS = vtkKWScubaLayer3DMRIS::New();
    layer3DMRIS->SetMRISProperties( this );

    layer = (vtkKWScubaLayer*)layer3DMRIS;

  }

  return layer;
}

void
vtkKWScubaLayerCollectionMRIS::LoadSurfaceFromFileName () {

  // Source object reads the surface and outputs a mesh.
  mSource = vtkFSSurfaceSource::New();
  mSource->MRISRead( mfnSurface.c_str() );
  mSource->Update();
}
