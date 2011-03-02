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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
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


#include <assert.h>
#include "vtkKWScubaLayerCollectionMRIS.h"
#include "vtkCommand.h"
#include "vtkDecimatePro.h"
#include "vtkFloatArray.h"
#include "vtkFSSurfaceSource.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkKWApplication.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWLabel.h"
#include "vtkKWLoadSaveDialog.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWLoadSaveButton.h"
#include "vtkKWProgressDialog.h"
#include "vtkKWScubaLayer2DMRIS.h"
#include "vtkKWScubaLayer3DMRIS.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtksys/SystemTools.hxx"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

extern "C" {
#include "colortab.h"
#include "mrisutils.h"
#include "error.h"
}

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollectionMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionMRIS, "$Revision: 1.4 $" );

vtkKWScubaLayerCollectionMRIS::vtkKWScubaLayerCollectionMRIS () :
  mfnSurface("") {

}

vtkKWScubaLayerCollectionMRIS::~vtkKWScubaLayerCollectionMRIS () {

}

void
vtkKWScubaLayerCollectionMRIS::SetSurfaceFileName ( const char* ifnSurface ) {

  // Set our file name and load the surface.
  mfnSurface = ifnSurface;
  this->LoadSurfaceFromFileName();
  assert( mSource.GetPointer() );

  // Find a good label based on the filename and set it in the
  // layer.
  string fnSurface = ifnSurface;
  string::size_type lastSlash = fnSurface.rfind( "/" );
  this->SetLabel( fnSurface.substr( lastSlash+1, string::npos ).c_str() );

  //
  // Decimator.
  //
  vtkSmartPointer<vtkKWProgressDialog> d = 
    vtkSmartPointer<vtkKWProgressDialog>::New();
  d->SetApplication( this->GetApplication() );
  d->SetWindowTitle( "Decimating Surface" );

  mDecimator = vtkSmartPointer<vtkDecimatePro>::New();
  mDecimator->AddObserver( vtkCommand::ProgressEvent, d );
  mDecimator->SetInputConnection( mSource->GetOutputPort() );
  mDecimator->SetTargetReduction( 0.9 );

  // Force the decimator to do its thing now.
  mDecimator->Update();
}

void
vtkKWScubaLayerCollectionMRIS::SetAnnotationAndColorTableFileNames 
( const char* ifnAnnotation, const char* ifnColorTable ) {
  
  // Set our file names and try loading the annotation.
  mfnAnnotationScalars = ifnAnnotation;
  if( ifnColorTable )
    mfnAnnotationTable = ifnColorTable;
  this->LoadAnnotationAndColorTableFromFileNames();

  if( mAnnotationScalars.GetPointer() && mAnnotationTable.GetPointer() ) {
    
    // Notify the layers that we now have scalar data.
    this->SendBroadcast( "DataAvailabilityChanged", NULL );
  }
}

void
vtkKWScubaLayerCollectionMRIS::AddControls ( vtkKWWidget* iPanel ) {

  vtkSmartPointer<vtkKWFrameWithLabel> frame = 
    vtkSmartPointer<vtkKWFrameWithLabel>::New();
  frame->SetParent( iPanel );
  frame->Create();
  frame->SetLabelText( "Annotation" );

  mLabelAnnotationFileName = vtkSmartPointer<vtkKWLabel>::New();
  mLabelAnnotationFileName->SetParent( frame->GetFrame() );
  mLabelAnnotationFileName->Create();
  mLabelAnnotationFileName->SetText( "No annotation loaded" );

  vtkSmartPointer<vtkKWPushButton> loadButton =
    vtkSmartPointer<vtkKWPushButton>::New();
  loadButton->SetParent( frame->GetFrame() );
  loadButton->Create();
  loadButton->SetText( "Load Annotation" );
  loadButton->SetCommand( this, "LoadAnnotationFromDlog" );
  
  this->Script( "pack %s %s -side top -fill x -expand y",
		mLabelAnnotationFileName->GetWidgetName(),
		loadButton->GetWidgetName() );

  this->Script( "pack %s -side top -fill x -expand y",
		frame->GetWidgetName() );
}

void
vtkKWScubaLayerCollectionMRIS::RemoveControls () {

  mLabelAnnotationFileName = NULL;
}

void
vtkKWScubaLayerCollectionMRIS::LoadAnnotationFromDlog () {

  // Create a Load dialog and set it up.
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load an annotation" );
  dialog->SetFileTypes( "{Annotation {.annot}} "
                        "{All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadAnnotation" );
  dialog->SetDefaultExtension( ".annot" );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadAnnotation" );

    string fnAnnotation = dialog->GetFileName();
    string fnColorTable = "";
    
    // We have an annotation file name, but we don't know if it has an
    // embedded color table. Check.
    int bHasColorTable = false;
    int err = 
      MRISisCTABPresentInAnnotation( fnAnnotation.c_str(), &bHasColorTable );
    if( err != ERROR_NONE ) {
      vtkKWMessageDialog::PopupMessage( this->GetApplication(), NULL,
					"Error", "Error loading the annotation file; possibly incorrect file format or incompatible version.",
					vtkKWMessageDialog::ErrorIcon );
    }
    
    // If it doesn't have one, ask for an LUT file.
    if( !bHasColorTable ) {
      
      dialog->SetTitle( "Choose a color lookup table" );
      dialog->SetFileTypes( "{LUT {.txt}} "
			    "{All {*}}" );
      dialog->RetrieveLastPathFromRegistry( "LoadLUT" );
      dialog->SetDefaultExtension( ".txt" );
      
      if( dialog->Invoke() ) {
	dialog->SaveLastPathToRegistry( "LoadLUT" );
	fnColorTable = dialog->GetFileName();
      }
    }
    
    // Set the file names. If we got a color table file name, pass it
    // in, otherwise pass NULL.
    this->SetAnnotationAndColorTableFileNames( fnAnnotation.c_str(),
		  fnColorTable != "" ? fnColorTable.c_str() : NULL );
  }
}

vtkFSSurfaceSource*
vtkKWScubaLayerCollectionMRIS::GetSource () const {
  return mSource.GetPointer();
}

vtkAlgorithmOutput*
vtkKWScubaLayerCollectionMRIS::GetNormalModeOutputPort () const {
  if( mSource.GetPointer() )
    return mSource->GetOutputPort();
  else 
    return NULL;
}

vtkAlgorithmOutput*
vtkKWScubaLayerCollectionMRIS::GetFastModeOutputPort () const {
  if( mDecimator.GetPointer() )
    return mDecimator->GetOutputPort();
  else 
    return NULL;
}

vtkPolyData*
vtkKWScubaLayerCollectionMRIS::GetNormalModeOutput () const {
  if( mSource.GetPointer() )
    return mSource->GetOutput();
  else 
    return NULL;
}

vtkPolyData*
vtkKWScubaLayerCollectionMRIS::GetFastModeOutput () const {
  if( mDecimator.GetPointer() )
    return mDecimator->GetOutput();
  else 
    return NULL;
}

vtkFloatArray*
vtkKWScubaLayerCollectionMRIS::GetScalarsValues () const {
  return mAnnotationScalars.GetPointer();
}

vtkScalarsToColors*
vtkKWScubaLayerCollectionMRIS::GetScalarsColors () const {
  return mAnnotationTable.GetPointer();
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

  try {

    // Source object reads the surface and outputs a mesh.
    mSource = vtkSmartPointer<vtkFSSurfaceSource>::New();
    mSource->MRISRead( mfnSurface.c_str() );
    mSource->Update();

  } catch ( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  }
}

void
vtkKWScubaLayerCollectionMRIS::LoadAnnotationAndColorTableFromFileNames () {

  int* lookupTable = NULL;
  COLOR_TABLE* table = NULL;
    
  try {
    
    if( NULL == mSource.GetPointer() )
      throw runtime_error( "Can't load annotation without first loading surface." );
    
    int cValues = mSource->GetNumberOfVertices();
    
    // Try to read the annotation.
    int eRead =
      MRISreadAnnotationIntoArray( mfnAnnotationScalars.c_str(), cValues,
				   &lookupTable );
    if( 0 != eRead )
      throw runtime_error ("Could not read annotation file");
    
    // Try to read a color table too if there is one.
    int bHasColorTable = false;
    eRead = MRISisCTABPresentInAnnotation( mfnAnnotationScalars.c_str(), 
					   &bHasColorTable );
    if( bHasColorTable ) {
      eRead =
	MRISreadCTABFromAnnotationIfPresent( mfnAnnotationScalars.c_str(), 
					     &table );
      if( ERROR_NONE != eRead ) {
	
	// If we didn't get a table, they have to have passed in a file
	// name for one.
	if( mfnAnnotationTable == "" ) {
	  throw runtime_error( "No color table present in overlay file; "
			       "need to specify an external file." );
	}
	
	//  Try to load the table.
	char* fnColors = strdup( mfnAnnotationTable.c_str() );
	table = CTABreadASCII( fnColors );
	free( fnColors );
	if( NULL == table )
	  throw runtime_error( string("Couldn't read color table ") + 
			       mfnAnnotationTable );
      }
    }

    // Init a float array.
    vtkSmartPointer<vtkFloatArray> scalars = 
      vtkSmartPointer<vtkFloatArray>::New();

    // Allocate our scalars.
    scalars->Allocate( cValues );
    scalars->SetNumberOfComponents( 1 );

    // Copy our array into the scalars.
    for( int nValue = 0; nValue < cValues; nValue ++ ) {
      int nEntry = 0;
      CTABfindAnnotation( table, lookupTable[nValue], &nEntry );
      scalars->InsertNextValue( static_cast<float>(nEntry) );
    }

    // Convert the CTAB to a vtkLookupTable.
    vtkSmartPointer<vtkFreesurferLookupTable> colors =
      vtkSmartPointer<vtkFreesurferLookupTable>::New();
    colors->BuildFromCTAB( table );

    // Save pointers.
    mAnnotationScalars = scalars;
    mAnnotationTable = colors;

    // Set the file name in the label.
    if( mLabelAnnotationFileName.GetPointer() ) 
      mLabelAnnotationFileName->SetText( vtksys::SystemTools::GetFilenameName( mfnAnnotationScalars.c_str() ).c_str() );

  }
  catch( exception& e ) {
    this->GetApplication()->ErrorMessage( e.what() );
  } 

  if( NULL != lookupTable ) free( lookupTable );
  if( NULL != table ) CTABfree( &table );
}
