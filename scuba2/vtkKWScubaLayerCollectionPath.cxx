/**
 * @file  vtkKWScubaLayerCollectionPath.cxx
 * @brief Implementation for path viewers.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.9 $
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

#include "vtkKWScubaLayerCollectionPath.h"
#include "vtkKWScubaLayer3DPath.h"
#include "vtkKWScubaLayer2DPath.h"
#include "vtkKWChangeColorButton.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButtonSetWithLabel.h"

#include "vtkFSVolumeSource.h"
#include "vtkObjectFactory.h"

#include "vtkSimplePointsReader.h"

#include "vtkImageReslice.h"
#include "vtkContourFilter.h"
#include "vtkTransform.h"

#include <string>
#include <iostream>
#include <vector>

using namespace std;

const double vtkKWScubaLayerCollectionPath::DEFAULT_COLOR[] = { 0.4, 0.5, 1.0 };

vtkStandardNewMacro( vtkKWScubaLayerCollectionPath );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionPath, "$Revision: 1.9 $" );

vtkKWScubaLayerCollectionPath::vtkKWScubaLayerCollectionPath ():
  mPathVolumeSource( NULL ),
  mSimplePointsReader( NULL ),
  mSamplesReader( NULL ),
  mfnPathVolume( "" ),
  mContourFilter( NULL ),
  mPathColorButton( NULL ) {

  mPathColor[ 0 ] = DEFAULT_COLOR[ 0 ];
  mPathColor[ 1 ] = DEFAULT_COLOR[ 1 ];
  mPathColor[ 2 ] = DEFAULT_COLOR[ 2 ];
  
  this->SetPathThresholdMode( PATH_THRESHOLD_LOW );
  
}

vtkKWScubaLayerCollectionPath::~vtkKWScubaLayerCollectionPath () {

  if( mContourFilter ) {
    mContourFilter->Delete();
  }
  
  if( mPathVolumeSource ) {
    mPathVolumeSource->Delete();
  }
  
  if( mSimplePointsReader ) {
    mSimplePointsReader->Delete();
  }
    
  if( mSamplesReader ) { 
    mSamplesReader->Delete();
  }
  
  if( !mInitialPaths.empty() ) {

    // this is a long for loop to iterate over the initial paths
    for( std::vector< std::vector< double* >* >::iterator itPaths = mInitialPaths.begin(); 
      itPaths != mInitialPaths.end(); itPaths++ ) {
        
      std::vector< double* >* path = *itPaths;
      
      // another long for loop to iterate over the points in a path
      for( std::vector< double* >::iterator itPoints = path->begin(); 
        itPoints != path->end(); itPoints++ ) {
          
        double *point = ( *itPoints );
        delete[] point;
        point = NULL;
        
      }
      
      path->clear();
      delete path;
      path = NULL;
      
    }
    
    mInitialPaths.clear();
        
  }
  
}

void 
vtkKWScubaLayerCollectionPath::SetVolumeFileName ( const char* ifnPathVolume ) {

  // Set our file name and load the volume.
  mfnPathVolume = ifnPathVolume;
  this->LoadPathVolumeFromFileName();

  // Find a good label based on the filename and set it in the
  // layer.
  string fnVolume = ifnPathVolume;
  string::size_type lastSlash = fnVolume.rfind( "/" );
  this->SetLabel( fnVolume.substr( lastSlash+1, string::npos ).c_str() );

}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionPath::GetPathVolumeSource () const {
  return mPathVolumeSource;
}

vtkSimplePointsReader* 
vtkKWScubaLayerCollectionPath::GetPathPointsSource () const {
  return mSimplePointsReader;
}

const std::vector< std::vector< double* >* > *
vtkKWScubaLayerCollectionPath::GetInitialPaths () const {
  return &mInitialPaths;
}


vtkPolyData* 
vtkKWScubaLayerCollectionPath::GetMesh () const {
  
  vtkPolyData* mesh = NULL;
  if( mContourFilter ) {
    mesh = mContourFilter->GetOutput();
  }
  
  return mesh;
};

void 
vtkKWScubaLayerCollectionPath::AddControls ( vtkKWWidget* iPanel ) {

  // create a color button for changing the color of the path
  mPathColorButton = vtkKWChangeColorButton::New();
  mPathColorButton->SetParent( iPanel );
  mPathColorButton->Create();
  
  // set the color to the default color
  mPathColorButton->SetColor( mPathColor[ 0 ], mPathColor[ 1 ], 
    mPathColor[ 2 ] );
  
  // set the callback
  mPathColorButton->SetCommand( this, "SetPathColor" );
  
  // Path Threshold
  vtkKWRadioButtonSetWithLabel* radBtnSetPathThreshold = vtkKWRadioButtonSetWithLabel::New();
  radBtnSetPathThreshold->SetParent( iPanel );
  radBtnSetPathThreshold->Create();
  radBtnSetPathThreshold->GetWidget()->PackHorizontallyOn();
  radBtnSetPathThreshold->SetLabelText( "Threshold: " );
  
  char sTclCmd[1024];

  vtkKWRadioButton* radBtnLowThreshold = radBtnSetPathThreshold->GetWidget()->AddWidget( 0 );
  radBtnLowThreshold->SetText( "Low" );
  sprintf( sTclCmd, "SetPathThresholdMode %d", PATH_THRESHOLD_LOW );
  radBtnLowThreshold->SetCommand( this, sTclCmd );
  if ( mPathThreshold == PATH_THRESHOLD_LOW )
    radBtnLowThreshold->SelectedStateOn();

  vtkKWRadioButton* radBtnHighThreshold = radBtnSetPathThreshold->GetWidget()->AddWidget( 1 );
  radBtnHighThreshold->SetText( "High" );
  sprintf( sTclCmd, "SetPathThresholdMode %d", PATH_THRESHOLD_HIGH );
  radBtnHighThreshold->SetCommand( this, sTclCmd );
  if ( mPathThreshold == PATH_THRESHOLD_HIGH )
    radBtnHighThreshold->SelectedStateOn();
    
  radBtnLowThreshold->Delete();
  radBtnHighThreshold->Delete();
    
  this->Script( "pack %s %s -side top -fill x",
                mPathColorButton->GetWidgetName(),
                radBtnSetPathThreshold->GetWidgetName() );
                  
}

void 
vtkKWScubaLayerCollectionPath::RemoveControls () {
  if( mPathColorButton ) {
    mPathColorButton->Delete();
    mPathColorButton = NULL;
  }
}

// called by the color change button too
void 
vtkKWScubaLayerCollectionPath::SetPathColor ( const double r, const double g, const double b ) {
  
  mPathColor[ 0 ] = r;
  mPathColor[ 1 ] = g;
  mPathColor[ 2 ] = b;

  this->SendBroadcast( "PathColorChanged", NULL );
    
}

void
vtkKWScubaLayerCollectionPath::GetPathColor( double &r, double &g, double &b ) const {
  r = mPathColor[ 0 ];
  g = mPathColor[ 1 ];
  b = mPathColor[ 2 ];
}

void 
vtkKWScubaLayerCollectionPath::GetPointColor( const int iPointIndex, double &r, double &g, double &b ) const {
  
  const int index = iPointIndex / 3;
  const int dimIndex = iPointIndex % 3;
  
  double* sample = mSamplesReader->GetOutput()->GetPoint( index );
    
  // let keep some all the time
  r = 0.6;
  g = sample[ dimIndex ];
  b = sample[ dimIndex ];
    
}

double
vtkKWScubaLayerCollectionPath::GetPointSampleValue( const int iPointIndex ) const {
  
  return mSamples[ iPointIndex ];
    
}

int
vtkKWScubaLayerCollectionPath::GetNumberOfSamples() const {
  return mSamples.size();
}

double 
vtkKWScubaLayerCollectionPath::GetPointProbabilityValue( const int iPointIndex ) const {
  return mProbabilities[ iPointIndex ];
}
    
int 
vtkKWScubaLayerCollectionPath::GetNumberOfProbabilities() const {
  return mProbabilities.size();
}

void 
vtkKWScubaLayerCollectionPath::SetPathThresholdMode ( int iMode ) {

  mPathThreshold = iMode;
  
  if( mContourFilter != NULL ) {
    mContourFilter->SetValue( 0, this->GetPathThreshold() );
    mContourFilter->Update();
  }
  
}

double 
vtkKWScubaLayerCollectionPath::GetPathThreshold () const {
  
  double threshold = 0.0;
  
  if( mPathThreshold == PATH_THRESHOLD_LOW ) {
    threshold = 0.00015;
  } else if ( mPathThreshold == PATH_THRESHOLD_HIGH ) {
    threshold = 0.0015;
  }
  
  return threshold;
  
}


vtkKWScubaLayer*
vtkKWScubaLayerCollectionPath::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {

  vtkKWScubaLayer* layer = NULL;

  if( vtkKWScubaView::TwoDee == iMode ) {

    vtkKWScubaLayer2DPath* layer2DPath = vtkKWScubaLayer2DPath::New();
    layer2DPath->SetPathProperties( this );
    
    layer = (vtkKWScubaLayer*)layer2DPath;
    
  } else if( vtkKWScubaView::ThreeDee == iMode ) {

    vtkKWScubaLayer3DPath* layer3DPath = vtkKWScubaLayer3DPath::New();
    layer3DPath->SetPathProperties( this );
    
    layer = (vtkKWScubaLayer*)layer3DPath;

  }

  return layer;
}

void 
vtkKWScubaLayerCollectionPath::LoadPathVolumeFromFileName () {

  // read in the coordinates of the path
  const char* sPathPoints = "/OptimalPath.txt";
  string fnOptimalPathPoints = this->GetFullFileName ( sPathPoints );
  mSimplePointsReader = vtkSimplePointsReader::New();
  mSimplePointsReader->SetFileName( fnOptimalPathPoints.c_str() );
  mSimplePointsReader->Update();
  
  // this will be a collection of paths -- only for DJ testing at the moment
  // NOTE from DJ: include this if you're testing out the initial paths
//  const char* sInitialPoints = "/InitialPath.txt";
//  string fnInitialPathPoints = this->GetFullFileName ( sInitialPoints );
//  this->ReadInitialPathPoints( fnInitialPathPoints.c_str() );
      
  // read in the values sampled along the path
  const char* sPathSamples = "/OptimalPathSamples.txt";
  string fnOptimalPathSamples = this->GetFullFileName ( sPathSamples );
  this->ReadScalars( fnOptimalPathSamples.c_str(), &mSamples );

  const char* sPathProbabilities = "/OptimalPathProbabilities.txt";  
  string fnOptimalPathProbabilities = this->GetFullFileName ( sPathProbabilities );
  this->ReadScalars( fnOptimalPathProbabilities.c_str(), &mProbabilities );
  
  // read in the path density
  mPathVolumeSource = vtkFSVolumeSource::New();
  mPathVolumeSource->MRIRead( mfnPathVolume.c_str() );
  mPathVolumeSource->Update();
  
  this->MakeMesh();
}

// Makes the path representation.
void 
vtkKWScubaLayerCollectionPath::MakeMesh () {
  //
  // Source object reads the volume and outputs structured points.
  //
  vtkFSVolumeSource* source = this->GetPathVolumeSource();
  
  // This transforms the voxel space source into RAS space.
  //
  vtkImageReslice* volumeToRAS = vtkImageReslice::New();
  volumeToRAS->SetInputConnection( source->GetOutputPort() );
  volumeToRAS->SetOutputDimensionality( 3 );

  // Get some values from the MRI.
  float RASBounds[6];
  source->GetRASBounds( RASBounds );
  
  // This rotates the volume to the proper orientation.
  double* rtv = source->GetRASToVoxelMatrix();
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->SetElement( 0, 0, rtv[0] );
  matrix->SetElement( 0, 1, rtv[1] );
  matrix->SetElement( 0, 2, rtv[2] );
  matrix->SetElement( 0, 3, 0 );
  matrix->SetElement( 1, 0, rtv[4] );
  matrix->SetElement( 1, 1, rtv[5] );
  matrix->SetElement( 1, 2, rtv[6] );
  matrix->SetElement( 1, 3, 0 );
  matrix->SetElement( 2, 0, rtv[8] );
  matrix->SetElement( 2, 1, rtv[9] );
  matrix->SetElement( 2, 2, rtv[10] );
  matrix->SetElement( 2, 3, 0 );
  matrix->SetElement( 3, 0, 0 );
  matrix->SetElement( 3, 1, 0 );
  matrix->SetElement( 3, 2, 0 );
  matrix->SetElement( 3, 3, 1 );
  
  vtkTransform* transform = vtkTransform::New();
  transform->SetMatrix( matrix );
  matrix->Delete();

  volumeToRAS->SetResliceTransform( transform );
  volumeToRAS->BorderOff();
  transform->Delete();

  // This sets our output extent.
  volumeToRAS->SetOutputExtent( (int)RASBounds[0], (int)RASBounds[1],
                                (int)RASBounds[2], (int)RASBounds[3],
                                (int)RASBounds[4], (int)RASBounds[5] );  
  
  mContourFilter = vtkContourFilter::New();
  mContourFilter->SetInputConnection( volumeToRAS->GetOutputPort() );
  mContourFilter->SetValue( 0, this->GetPathThreshold() );
  
  volumeToRAS->Delete();
    
}

string
vtkKWScubaLayerCollectionPath::GetFullFileName ( const char* sShortFileName ) const {
  
  // part of the string to be replaced
  const char* sReplace = "/OptimalPathDensity";

  string fullFileName = mfnPathVolume;
  size_t nFound = fullFileName.find( sReplace );

  if( nFound != string::npos ) {
    const int amountToReplace = fullFileName.size() - nFound;
    fullFileName.replace( nFound, amountToReplace, sShortFileName );
  }
  
  return fullFileName;
  
}

void 
vtkKWScubaLayerCollectionPath::ReadScalars( const char* fnScalar, 
  std::vector< double > *ioScalars ) {

  string s;
  ifstream scalarsFile( fnScalar );

  double scalar = 0.0;
  
  // read all the samples
  while( scalarsFile >> s ) {
    
    // convert the sample to double format
    stringstream stream;
    stream << s;
    stream >> scalar;
    
    // add the sample the samples vector
    ioScalars->insert( ioScalars->end(), scalar );
  }
  
}

void 
vtkKWScubaLayerCollectionPath::ReadInitialPathPoints( const char* fnInitialPathPoints ) {
    
  // this is the delimiter for indicating the start of a new pathway
  const double pathDelimiter = -999.0;
  
  vtkSimplePointsReader *initialPointsReader = vtkSimplePointsReader::New();
  initialPointsReader->SetFileName( fnInitialPathPoints );
  initialPointsReader->Update();
    
  const int cTotalPoints = initialPointsReader->GetOutput()->GetNumberOfPoints();
    
  bool bCreateNewPath = true;
    
  std::vector< double* > *path = NULL;
  
  // go through all the points and at create new tubes when needed
  for( int nTotalPoints=0; nTotalPoints<cTotalPoints; nTotalPoints++ ) {

    const double *point = initialPointsReader->GetOutput()->GetPoint( nTotalPoints );
    
    // determine if this is a new path
     if( !bCreateNewPath ) {
      
      bCreateNewPath = true;
      
      // if all three columns are -999, then we've a new tube
      for( int nDim = 0; nDim < 3 && bCreateNewPath; nDim++ ) {
        
        if( pathDelimiter != point[ nDim ] ) {
          bCreateNewPath = false;
        }
        
      }
      
      // if we've made it through the above loop and all three elements are the
      // new tube delimiter, then we've got a new tube
      
    }
    
    // create the proper storage if we need to
    if( bCreateNewPath ) {
      // allocate memory for a new path
      path = new std::vector< double* >();
      
      // push the path
      mInitialPaths.push_back( path );
      
      bCreateNewPath = false;
    } else {
      // we don't want to copy the delimiter over, so we skip this if we need
      // to create a new path
          
      // copy the point and save it to the path
      double *pathPoints = new double[ 3 ];
      for( int nDim=0; nDim<3; nDim++ ) {
        pathPoints[ nDim ] = point[ nDim ];
      }
      path->push_back( pathPoints );
    }
    
  }  
  
  initialPointsReader->Delete();
  
}
