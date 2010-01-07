/**
 * @file  LayerMRI.cpp
 * @brief Layer class for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/01/07 23:33:04 $
 *    $Revision: 1.49 $
 *
 * Copyright (C) 2008-2009,
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

#include <wx/wx.h>
#include "LayerMRI.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkCubeSource.h"
#include "vtkVolume.h"
#include "vtkImageThreshold.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkContourFilter.h"
#include "vtkTubeFilter.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkAppendPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

#include "LayerPropertiesMRI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "MainWindow.h"
#include "vtkMath.h"
#include "BuildContourThread.h"


#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

LayerMRI::LayerMRI( LayerMRI* ref ) : LayerVolumeBase(),
    m_volumeSource( NULL),
    m_volumeRef( ref ? ref->GetSourceVolume() : NULL ),
    m_bResampleToRAS( false ),
    m_bReorient( false ),
    m_nSampleMethod( SAMPLE_NEAREST )
{
  m_strTypeNames.push_back( "MRI" );

  for ( int i = 0; i < 3; i++ )
  {
    // m_nSliceNumber[i] = 0;
    m_sliceActor2D[i] = vtkImageActor::New();
    m_sliceActor3D[i] = vtkImageActor::New();
    /* 
    m_sliceActor2D[i]->GetProperty()->SetAmbient( 1 );
    m_sliceActor2D[i]->GetProperty()->SetDiffuse( 0 );
    m_sliceActor3D[i]->GetProperty()->SetAmbient( 1 );
    m_sliceActor3D[i]->GetProperty()->SetDiffuse( 0 );
    */
    m_sliceActor2D[i]->InterpolateOff();
    m_sliceActor3D[i]->InterpolateOff();
    
    m_glyphActor2D[i] = vtkActor::New();
    m_glyphActor3D[i] = vtkActor::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_glyphActor2D[i]->SetMapper( mapper );
    m_glyphActor3D[i]->SetMapper( mapper2 );
  }
  mProperties = new LayerPropertiesMRI();
  mProperties->AddListener( this );

  m_actorContour = vtkActor::New();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
  m_actorContour->SetMapper( mapper );
  m_propVolume = vtkVolume::New();
  
  m_nThreadID = 0;
  
  private_buf1_3x3 = new double*[3];
  private_buf2_3x3 = new double*[3];
  for ( int i = 0; i < 3; i++ )
  {
    private_buf1_3x3[i] = new double[3];
    private_buf2_3x3[i] = new double[3];
  }
}

LayerMRI::~LayerMRI()
{	
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->Delete();
		m_sliceActor3D[i]->Delete();
    m_glyphActor2D[i]->Delete();
    m_glyphActor3D[i]->Delete();
	}
	m_actorContour->Delete();
	m_propVolume->Delete();
	
  if ( m_sFilename.size() > 0 )
    mProperties->SaveSettings( m_sFilename.c_str() );
	delete mProperties;
	
	if ( m_volumeSource )
		delete m_volumeSource;
    
  for ( size_t i = 0; i < m_segActors.size(); i++ )
  {
    m_segActors[i].actor->Delete();
  }
  
  for ( int i = 0; i < 3; i++ )
  {
    delete[] private_buf1_3x3[i];
    delete[] private_buf2_3x3[i];
  }
  delete[] private_buf1_3x3;
  delete[] private_buf2_3x3;
}

void LayerMRI::SetResampleToRAS( bool bResample )
{
  m_bResampleToRAS = bResample;
}

bool LayerMRI::LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_volumeSource )
    delete m_volumeSource;

  m_volumeSource = new FSVolume( m_volumeRef );
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_volumeSource->SetInterpolationMethod( m_nSampleMethod );

  if ( !m_volumeSource->MRIRead(  m_sFilename.c_str(),
                                  m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL,
                                  wnd,
                                  event ) )
    return false;

  event.SetInt( 100 );
  wxPostEvent( wnd, event );

  InitializeVolume();
  InitializeActors();

  mProperties->SetVolumeSource( m_volumeSource );
  mProperties->RestoreSettings( m_sFilename.c_str() );
  
  return true;
}

bool LayerMRI::Create( LayerMRI* mri, bool bCopyVoxelData, int data_type )
{
  if ( m_volumeSource )
    delete m_volumeSource;

  m_volumeSource = new FSVolume( mri->m_volumeSource );
  if ( !m_volumeSource->Create( mri->m_volumeSource, bCopyVoxelData, data_type ) )
    return false;

  m_bResampleToRAS = mri->m_bResampleToRAS;
  m_bReorient = mri->m_bReorient;
  m_imageDataRef = mri->GetImageData();
// if ( m_imageDataRef.GetPointer() )
  if ( m_imageDataRef != NULL )
  {
    SetWorldOrigin( mri->GetWorldOrigin() );
    SetWorldVoxelSize( mri->GetWorldVoxelSize() );
    SetWorldSize( mri->GetWorldSize() );

    m_imageData = m_volumeSource->GetImageOutput();

    InitializeActors();

    mProperties->SetVolumeSource( m_volumeSource );
    // mProperties->SetColorMap( LayerPropertiesMRI::LUT );

    m_sFilename = "";

    if ( bCopyVoxelData )
    {
      mProperties->CopySettings( mri->mProperties ); 
      SetModified();
    }
  }

  return true;
}

void LayerMRI::SetReorient( bool bReorient )
{
  m_bReorient = bReorient;
}

bool LayerMRI::SaveVolume( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_sFilename.size() == 0 || m_imageData == NULL )
    return false;
  
  /*
  if ( IsModified() || m_bReorient )
    m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event, !m_bReorient );
  else
  {
    if ( !m_volumeSource->Restore( m_sFilename.c_str(), 
                                  m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL, 
                                  wnd, 
                                  event ) )
    {
      cerr << "Restore from converted volume." << endl;
      m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event, !m_bReorient );
    }
  }
  */
  if ( !m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event, !m_bReorient ) )
    return false;
  
  // now first save a copy of the old file if required
  if ( MainWindow::GetMainWindowPointer()->GetSaveCopy() )
  {
    wxString new_fn = m_sFilename + "~";
    rename( m_sFilename.c_str(), new_fn.c_str() );
  }
  
  bool bSaved = m_volumeSource->MRIWrite( m_sFilename.c_str(), !m_bReorient );
  m_bModified = !bSaved;

  event.SetInt( 99 );
  wxPostEvent( wnd, event );

  return bSaved;
}

bool LayerMRI::Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  m_bResampleToRAS = false;
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  if ( IsModified() )
  {
    if ( !m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event ) )
      return false;
  }
  else
  {
    if ( !m_volumeSource->Restore( m_sFilename.c_str(), 
          m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL, 
          wnd, 
          event ) )
    {
      cerr << "Failed to load original volume. Restore from converted volume (not good)." << endl;
      m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event, !m_bReorient );
    }
  }

  int nSampleMethod = rotations[0].SampleMethod;
  if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
    nSampleMethod = SAMPLE_NEAREST;
  bool ret = m_volumeSource->Rotate( rotations, wnd, event, nSampleMethod );
  m_imageData = m_volumeSource->GetImageOutput();
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->SetInput( m_imageData );
  }
  
  return ret;
}

void LayerMRI::InitializeVolume()
{
  if ( m_volumeSource == NULL )
    return;

  FSVolume* source = m_volumeSource;

  float RASBounds[6];
  source->GetBounds( RASBounds );
  m_dWorldOrigin[0] = RASBounds[0];
  m_dWorldOrigin[1] = RASBounds[2];
  m_dWorldOrigin[2] = RASBounds[4];
  source->GetPixelSize( m_dWorldVoxelSize );

  m_dWorldSize[0] = ( ( int )( (RASBounds[1] - RASBounds[0]) / m_dWorldVoxelSize[0] ) ) * m_dWorldVoxelSize[0];
  m_dWorldSize[1] = ( ( int )( (RASBounds[3] - RASBounds[2]) / m_dWorldVoxelSize[1] ) ) * m_dWorldVoxelSize[1];
  m_dWorldSize[2] = ( ( int )( (RASBounds[5] - RASBounds[4]) / m_dWorldVoxelSize[2] ) ) * m_dWorldVoxelSize[2];

  m_imageData = source->GetImageOutput();
}


void LayerMRI::InitializeActors()
{
  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
    mReslice[i]->SetInput( m_imageData );
//  mReslice[i]->SetOutputSpacing( sizeX, sizeY, sizeZ );
    mReslice[i]->BorderOff();

    // This sets us to extract slices.
    mReslice[i]->SetOutputDimensionality( 2 );

    // This will change depending what orienation we're in.
    mReslice[i]->SetResliceAxesDirectionCosines( 1, 0, 0,
        0, 1, 0,
        0, 0, 1 );

    // This will change to select a different slice.
    mReslice[i]->SetResliceAxesOrigin( 0, 0, 0 );

    //
    // Image to colors using color table.
    //
    mColorMap[i] = vtkSmartPointer<vtkImageMapToColors>::New();
    mColorMap[i]->SetLookupTable( mProperties->GetGrayScaleTable() );
    mColorMap[i]->SetInput( mReslice[i]->GetOutput() );
    mColorMap[i]->SetOutputFormatToRGBA();
    mColorMap[i]->PassAlphaToOutputOn();

    //
    // Prop in scene with plane mesh and texture.
    //
    m_sliceActor2D[i]->SetInput( mColorMap[i]->GetOutput() );
    m_sliceActor3D[i]->SetInput( mColorMap[i]->GetOutput() );

    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }

  this->UpdateResliceInterpolation();
  this->UpdateTextureSmoothing();
  this->UpdateOpacity();
  this->UpdateColorMap();
}

void LayerMRI::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetOpacity( mProperties->GetOpacity() );
    m_sliceActor3D[i]->SetOpacity( mProperties->GetOpacity() );
  }
}

void LayerMRI::UpdateColorMap ()
{
  assert( mProperties );

  for ( int i = 0; i < 3; i++ )
    mColorMap[i]->SetActiveComponent( m_nActiveFrame );

  switch ( mProperties->GetColorMap() )
  {
  case LayerPropertiesMRI::NoColorMap:
    for ( int i = 0; i < 3; i++ )
      mColorMap[i]->SetLookupTable( NULL );
    break;
  case LayerPropertiesMRI::Grayscale:
    for ( int i = 0; i < 3; i++ )
      mColorMap[i]->SetLookupTable( mProperties->GetGrayScaleTable() );
    m_actorContour->GetMapper()->SetLookupTable( mProperties->GetGrayScaleTable() );
    break;
  case LayerPropertiesMRI::Heat:
    for ( int i = 0; i < 3; i++ )
      mColorMap[i]->SetLookupTable( mProperties->GetHeatScaleTable() );
    m_actorContour->GetMapper()->SetLookupTable( mProperties->GetHeatScaleTable() );
    break;
  case LayerPropertiesMRI::Jet:
    for ( int i = 0; i < 3; i++ )
      mColorMap[i]->SetLookupTable( mProperties->GetJetScaleTable() );
    m_actorContour->GetMapper()->SetLookupTable( mProperties->GetJetScaleTable() );
    break;
  case LayerPropertiesMRI::LUT:
    for ( int i = 0; i < 3; i++ )
      mColorMap[i]->SetLookupTable( mProperties->GetLUTTable() );
    m_actorContour->GetMapper()->SetLookupTable( mProperties->GetLUTTable() );
    break;
  default:
    break;
  }
}

void LayerMRI::UpdateResliceInterpolation ()
{
  assert( mProperties );

  for ( int i = 0; i < 3; i++ )
  {
    if ( mReslice[i].GetPointer() )
    {
      mReslice[i]->SetInterpolationMode( mProperties->GetResliceInterpolation() );
    }
  }
}

void LayerMRI::UpdateTextureSmoothing ()
{
  assert( mProperties );

  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetInterpolate( mProperties->GetTextureSmoothing() );
    m_sliceActor3D[i]->SetInterpolate( mProperties->GetTextureSmoothing() );
  }
}

void LayerMRI::UpdateContour( int nSegValue )
{  
  UpdateContourActor( nSegValue );
}

void LayerMRI::UpdateContourActor( int nSegValue )
{
  /*
  double dTh1 = GetProperties()->GetContourMinThreshold();
  double dTh2 = GetProperties()->GetContourMaxThreshold();
  if ( nSegValue >= 0 )
  {
    dTh1 = nSegValue - 0.5;
    dTh2 = nSegValue + 0.5;
  }
  MyUtils::BuildContourActor( GetImageData(), dTh1, dTh2, m_actorContour );
  */
  
  // Generate a new thread id before creating the thread. so that mainwindow will be able to determine 
  // if a build contour result is already expired, by comparing the returned id and current id. If they
  // are different, it means a new thread is rebuilding the contour
  m_nThreadID++;
  BuildContourThread* thread = new BuildContourThread( MainWindow::GetMainWindowPointer() );
  thread->BuildContour( this, nSegValue, m_nThreadID );
}

// Contour mapper is ready, attach it to the actor   
void LayerMRI::RealizeContourActor()
{
  if ( m_actorContourTemp.GetPointer() && m_actorContourTemp->GetMapper() )
  {
    m_actorContour->SetMapper( m_actorContourTemp->GetMapper() );
    UpdateContourColor();
  }
}

void LayerMRI::ShowContour()
{
  if ( !m_actorContourTemp.GetPointer() )
    UpdateContour();
}

void LayerMRI::UpdateContourColor()
{
  if ( mProperties->GetContourUseImageColorMap() )
  {
    m_actorContour->GetMapper()->ScalarVisibilityOn();
  }
  else
  {
    m_actorContour->GetMapper()->ScalarVisibilityOff();
    m_actorContour->GetProperty()->SetColor( mProperties->GetContourColor() );
  }
  UpdateColorMap();
}

void LayerMRI::UpdateVolumeRendering()
{
  /*
  if ( GetProperties()->GetShowAsContour() )
  {
    MyUtils::BuildVolume( GetImageData(),
                          GetProperties()->GetContourMinThreshold(),
                          GetProperties()->GetContourMaxThreshold(),
                          m_propVolume );
  }
  */
}

void LayerMRI::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  wxASSERT ( nPlane >= 0 && nPlane <= 2 );

  if ( mProperties->GetDisplayVector() || mProperties->GetDisplayTensor() )
    renderer->AddViewProp( m_glyphActor2D[nPlane] );
  else
    renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerMRI::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
	bool bContour = GetProperties()->GetShowAsContour();
	for ( int i = 0; i < 3; i++ )
	{
		if ( !bContour && ( bSliceVisibility == NULL || bSliceVisibility[i] ) )
    {
      if ( mProperties->GetDisplayVector() || mProperties->GetDisplayTensor() )
        renderer->AddViewProp( m_glyphActor3D[i] );
			else
        renderer->AddViewProp( m_sliceActor3D[i] );
    } 
	}
	
	if ( bContour )
	{
   // for ( size_t i = 0; i < m_segActors.size(); i ++ )
	 // renderer->AddViewProp( m_segActors[i].actor );
    renderer->AddViewProp( m_actorContour );
	}
}

/*
void LayerMRI::SetSliceNumber( int* sliceNumber )
{
 if ( sliceNumber[0] != m_nSliceNumber[0] || sliceNumber[1] != m_nSliceNumber[1] ||
   sliceNumber[2] != m_nSliceNumber[2] )
 {
  m_nSliceNumber[0] = sliceNumber[0];
  m_nSliceNumber[1] = sliceNumber[1];
  m_nSliceNumber[2] = sliceNumber[2];


 }
}
*/

void LayerMRI::SetSlicePositionToWorldCenter()
{
  if ( m_volumeSource == NULL )
    return;

  // Get some values from the MRI.
  double pos[3];
  for ( int i = 0; i < 3; i++ )
    pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];

  SetSlicePosition( pos );
}

void LayerMRI::OnSlicePositionChanged( int nPlane )
{
  if ( m_volumeSource == NULL || nPlane < 0 || nPlane > 2)
    return;

  assert( mProperties );
 
  // display slice image
  vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  switch ( nPlane )
  {
    case 0:
      m_sliceActor2D[0]->PokeMatrix( matrix );
      m_sliceActor2D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
      m_sliceActor2D[0]->RotateX( 90 );
      m_sliceActor2D[0]->RotateY( -90 );
      m_sliceActor3D[0]->PokeMatrix( matrix );
      m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
      m_sliceActor3D[0]->RotateX( 90 );
      m_sliceActor3D[0]->RotateY( -90 );
  
      // Putting negatives in the reslice axes cosines will flip the
      // image on that axis.
      mReslice[0]->SetResliceAxesDirectionCosines( 0, -1, 0,
          0, 0, 1,
          1, 0, 0 );
      mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
      mReslice[0]->Modified();
      break;
    case 1:
      m_sliceActor2D[1]->PokeMatrix( matrix );
      m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
      m_sliceActor2D[1]->RotateX( 90 );
      // m_sliceActor2D[1]->RotateY( 180 );
      m_sliceActor3D[1]->PokeMatrix( matrix );
      m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
      m_sliceActor3D[1]->RotateX( 90 );
      // m_sliceActor3D[1]->RotateY( 180 );
  
      // Putting negatives in the reslice axes cosines will flip the
      // image on that axis.
      mReslice[1]->SetResliceAxesDirectionCosines( 1, 0, 0,
          0, 0, 1,
          0, 1, 0 );
      mReslice[1]->SetResliceAxesOrigin( 0, m_dSlicePosition[1], 0 );
      mReslice[1]->Modified();
      break;
    case 2:
      m_sliceActor2D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
      // m_sliceActor2D[2]->RotateY( 180 );
      m_sliceActor3D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
      // m_sliceActor3D[2]->RotateY( 180 );
  
      mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
          0, 1, 0,
          0, 0, 1 );
      mReslice[2]->SetResliceAxesOrigin( 0, 0, m_dSlicePosition[2]  );
      mReslice[2]->Modified();
      break;
  }
  // display 4D data as vector
  if ( mProperties->GetDisplayVector() )
  {
    UpdateVectorActor( nPlane );
  }
  else if ( mProperties->GetDisplayTensor() )
  {
    UpdateTensorActor( nPlane );
  }
}

LayerPropertiesMRI* LayerMRI::GetProperties()
{
  return mProperties;
}

void LayerMRI::DoListenToMessage( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->UpdateColorMap();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "ResliceInterpolationChanged" )
  {
    this->UpdateResliceInterpolation();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "OpacityChanged" )
  {
    this->UpdateOpacity();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "TextureSmoothingChanged" )
  {
    this->UpdateTextureSmoothing();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "LayerContourChanged" )
  {
    if ( this->GetProperties()->GetShowAsContour() )
      this->UpdateContour();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "LayerContourColorChanged" )
  {
    this->UpdateContourColor();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "LayerContourShown" )
  {
    if ( this->GetProperties()->GetShowAsContour() )
      this->ShowContour();
    this->SendBroadcast( iMessage, this, this );
  }
  else if ( iMessage == "DisplayModeChanged" )
  {
    for ( int i = 0; i < 3; i++ )
    {
      if ( mProperties->GetDisplayTensor() && 
           mProperties->GetTensorRepresentation() == LayerPropertiesMRI::TR_Ellipsoid )
      {
        m_glyphActor2D[i]->GetProperty()->SetInterpolationToGouraud();
        m_glyphActor3D[i]->GetProperty()->SetInterpolationToGouraud();
      }
      else
      {
        m_glyphActor2D[i]->GetProperty()->SetInterpolationToFlat();
        m_glyphActor3D[i]->GetProperty()->SetInterpolationToFlat();
      }
    }
    if ( mProperties->GetDisplayVector() )
      this->UpdateVectorActor();
    else if ( mProperties->GetDisplayTensor() )
      this->UpdateTensorActor();
    this->SendBroadcast( iMessage, this, this );
  }
}

void LayerMRI::SetVisible( bool bVisible )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_glyphActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_glyphActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
  }
  m_actorContour->SetVisibility( bVisible ? 1 : 0 );
  this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerMRI::IsVisible()
{
  return m_sliceActor2D[0]->GetVisibility() > 0;
}

double LayerMRI::GetVoxelValue( double* pos )
{ 
  if ( m_imageData == NULL )
    return 0;

  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  int* ext = m_imageData->GetExtent();

  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos[i] - orig[i] ) / vsize[i] + 0.5 );
  }

  if ( n[0] < ext[0] || n[0] > ext[1] ||
       n[1] < ext[2] || n[1] > ext[3] ||
       n[2] < ext[4] || n[2] > ext[5] )
    return 0;
  else
    return m_imageData->GetScalarComponentAsDouble( n[0], n[1], n[2], m_nActiveFrame );
}

double LayerMRI::GetVoxelValueByOriginalIndex( int i, int j, int k )
{
  return m_volumeSource->GetVoxelValue( i, j, k, m_nActiveFrame );
}

void LayerMRI::SetModified()
{
  mReslice[0]->Modified();
  mReslice[1]->Modified();
  mReslice[2]->Modified();

  LayerVolumeBase::SetModified();
}

std::string LayerMRI::GetLabelName( double value )
{
  int nIndex = (int)value;
  if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
  {
    COLOR_TABLE* ct = GetProperties()->GetLUTCTAB();
    if ( !ct )
      return "";
    char name[128];
    int nValid = 0;
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    if ( nIndex > 0 && nIndex < nTotalCount )
    {
      CTABisEntryValid( ct, nIndex, &nValid );
      if ( nValid && CTABcopyName( ct, nIndex, name, 128 ) == 0 )
      {
        return name;
      }
    }
  }

  return "";
}

void LayerMRI::RemapPositionToRealRAS( const double* pos_in, double* pos_out )
{
  m_volumeSource->TargetToRAS( pos_in, pos_out );
}

void LayerMRI::RemapPositionToRealRAS( double x_in, double y_in, double z_in,
                                       double& x_out, double& y_out, double& z_out )
{
  m_volumeSource->TargetToRAS( x_in, y_in, z_in, x_out, y_out, z_out );
}

void LayerMRI::RASToTarget( const double* pos_in, double* pos_out )
{
  m_volumeSource->RASToTarget( pos_in, pos_out );
}

int LayerMRI::GetNumberOfFrames()
{
  if ( m_imageData )
    return m_imageData->GetNumberOfScalarComponents();
  else
    return 1;
}

void LayerMRI::SetActiveFrame( int nFrame )
{
  if ( nFrame != m_nActiveFrame && nFrame >= 0 && nFrame < this->GetNumberOfFrames() )
  {
    m_nActiveFrame = nFrame;
    this->DoListenToMessage( "ColorMapChanged", this, this );
    this->SendBroadcast( "LayerActiveFrameChanged", this, this );
  }
}

bool LayerMRI::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( m_sliceActor3D[i] == prop )
      return true;
  }
  return false;
}

void LayerMRI::RASToOriginalIndex( const double* pos, int* n )
{
  m_volumeSource->RASToOriginalIndex( (float)(pos[0]), (float)(pos[1]), (float)(pos[2]),
                                      n[0], n[1], n[2] );
}

void LayerMRI::OriginalIndexToRAS( const int* n, double* pos )
{
  float x, y, z;
  m_volumeSource->OriginalIndexToRAS( n[0], n[1], n[2], x, y, z );
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
}

void LayerMRI::UpdateVectorActor()
{
  for ( int i = 0; i < 3; i++ )
  {
    UpdateVectorActor( i );
  }
}

void LayerMRI::UpdateVectorActor( int nPlane )
{
  UpdateVectorActor( nPlane, m_imageData );
}

void LayerMRI::UpdateVectorActor( int nPlane, vtkImageData* imagedata )
{
  double* pos = GetSlicePosition();
  double* orig = imagedata->GetOrigin();
  double* voxel_size = imagedata->GetSpacing();
  int* dim = imagedata->GetDimensions();
  int n[3];
  for ( int i = 0; i < 3; i++ )
    n[i] = (int)( ( pos[i] - orig[i] ) / voxel_size[i] + 0.5 );

//  vtkPolyData* polydata = vtkPolyDataMapper::SafeDownCast( m_vectorActor2D[nPlane]->GetMapper() )->GetInput();
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  scalars->SetNumberOfComponents( 4 );  
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  polydata->GetPointData()->SetScalars( scalars );
  if ( n[0] < 0 || n[0] >= dim[0] ||
       n[1] < 0 || n[1] >= dim[1] ||
       n[2] < 0 || n[2] >= dim[2] )
  {
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
    return;
  }
  
  int nCnt = 0;
  double scale = min( min( voxel_size[0], voxel_size[1] ), voxel_size[2] ) / 1.8;
  
  if ( GetProperties()->GetVectorRepresentation() == LayerPropertiesMRI::VR_Bar )
  {
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInput( polydata );
    tube->SetNumberOfSides( 4 );
    tube->SetRadius( scale * 0.3 );
    tube->CappingOn();
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( tube->GetOutput() );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( tube->GetOutput() );
  }
  else
  {
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
  }
  
  unsigned char c[4] = { 0, 0, 0, 255 };
  switch ( nPlane )
  {
    case 0:
      for ( int i = 0; i < dim[1]; i++ )
      {
        for ( int j = 0; j < dim[2]; j++ )
        {
          double v[3];
          double pt[3];
          v[0] = imagedata->GetScalarComponentAsDouble( n[0], i, j, 0 );
          v[1] = imagedata->GetScalarComponentAsDouble( n[0], i, j, 1 );
          v[2] = imagedata->GetScalarComponentAsDouble( n[0], i, j, 2 );
          if ( vtkMath::Normalize( v ) != 0 )
          {
            v[1] = -v[1];         // by default invert Y !!
            if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_X )
              v[0] = -v[0];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Y )
              v[1] = -v[1];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Z )
              v[2] = -v[2];
        
            pt[0] = orig[0] + voxel_size[0] * n[0];
            pt[1] = orig[1] + voxel_size[1] * i;
            pt[2] = orig[2] + voxel_size[2] * j;
            lines->InsertNextCell( 2 );
            points->InsertNextPoint( pt[0] + scale * v[0], 
                                    pt[1] + scale * v[1], 
                                    pt[2] + scale * v[2] );
            points->InsertNextPoint( pt[0] - scale * v[0], 
                                    pt[1] - scale * v[1], 
                                    pt[2] - scale * v[2] );
            lines->InsertCellPoint( nCnt++ );
            lines->InsertCellPoint( nCnt++ );
            c[0] = (int)(fabs( v[0] *255 ) );
            c[1] = (int)(fabs( v[1] *255 ) );
            c[2] = (int)(fabs( v[2] *255 ) );
            scalars->InsertNextTupleValue( c );
            scalars->InsertNextTupleValue( c );
          }
        }
      }
      break;
    case 1:
      for ( int i = 0; i < dim[0]; i++ )
      {
        for ( int j = 0; j < dim[2]; j++ )
        {
          double v[3];
          double pt[3];
          v[0] = imagedata->GetScalarComponentAsDouble( i, n[1], j, 0 );
          v[1] = imagedata->GetScalarComponentAsDouble( i, n[1], j, 1 );
          v[2] = imagedata->GetScalarComponentAsDouble( i, n[1], j, 2 );
          if ( vtkMath::Normalize( v ) != 0 )
          {
            v[1] = -v[1];         // by default invert Y !!
            if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_X )   
              v[0] = -v[0];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Y )
              v[1] = -v[1];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Z )
              v[2] = -v[2];
            
            pt[0] = orig[0] + voxel_size[0] * i;
            pt[1] = orig[1] + voxel_size[1] * n[1];
            pt[2] = orig[2] + voxel_size[2] * j;
            lines->InsertNextCell( 2 );
            points->InsertNextPoint( pt[0] + scale * v[0], 
                                    pt[1] + scale * v[1], 
                                    pt[2] + scale * v[2] );
            points->InsertNextPoint( pt[0] - scale * v[0], 
                                    pt[1] - scale * v[1], 
                                    pt[2] - scale * v[2] );
            lines->InsertCellPoint( nCnt++ );
            lines->InsertCellPoint( nCnt++ );
            c[0] = (int)(fabs( v[0] *255 ) );
            c[1] = (int)(fabs( v[1] *255 ) );
            c[2] = (int)(fabs( v[2] *255 ) );
            scalars->InsertNextTupleValue( c );
            scalars->InsertNextTupleValue( c );
          }
        }
      }
      break;
    case 2:
      for ( int i = 0; i < dim[0]; i++ )
      {
        for ( int j = 0; j < dim[1]; j++ )
        {
          double v[3];
          double pt[3];
          v[0] = imagedata->GetScalarComponentAsDouble( i, j, n[2], 0 );
          v[1] = imagedata->GetScalarComponentAsDouble( i, j, n[2], 1 );
          v[2] = imagedata->GetScalarComponentAsDouble( i, j, n[2], 2 );
          if ( vtkMath::Normalize( v ) != 0 )
          {
            v[1] = -v[1];         // by default invert Y !!
            if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_X )   
              v[0] = -v[0];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Y )
              v[1] = -v[1];
            else if ( GetProperties()->GetVectorInversion() == LayerPropertiesMRI::VI_Z )
              v[2] = -v[2];
            
            pt[0] = orig[0] + voxel_size[0] * i;
            pt[1] = orig[1] + voxel_size[1] * j;
            pt[2] = orig[2] + voxel_size[2] * n[2];
            lines->InsertNextCell( 2 );
            points->InsertNextPoint( pt[0] + scale * v[0], 
                                    pt[1] + scale * v[1], 
                                    pt[2] + scale * v[2] );
            points->InsertNextPoint( pt[0] - scale * v[0], 
                                    pt[1] - scale * v[1], 
                                    pt[2] - scale * v[2] );
            lines->InsertCellPoint( nCnt++ );
            lines->InsertCellPoint( nCnt++ );
            c[0] = (int)(fabs( v[0] *255 ) );
            c[1] = (int)(fabs( v[1] *255 ) );
            c[2] = (int)(fabs( v[2] *255 ) );
            scalars->InsertNextTupleValue( c );
            scalars->InsertNextTupleValue( c );
          }
        }
      }
      break;
    default:
      break;
  }
}

void LayerMRI::UpdateTensorActor()
{
  for ( int i = 0; i < 3; i++ )
  {
    UpdateTensorActor( i );
  }
}

void LayerMRI::UpdateTensorActor( int nPlane, vtkImageData* imagedata_in )
{
  vtkImageData* imagedata = imagedata_in;
  if ( !imagedata )
    imagedata = m_imageData;
  
  double* pos = GetSlicePosition();
  double* orig = imagedata->GetOrigin();
  double* voxel_size = imagedata->GetSpacing();
  int* dim = imagedata->GetDimensions();
  int n[3];
  for ( int i = 0; i < 3; i++ )
    n[i] = (int)( ( pos[i] - orig[i] ) / voxel_size[i] + 0.5 );

  if ( n[0] < 0 || n[0] >= dim[0] ||
       n[1] < 0 || n[1] >= dim[1] ||
       n[2] < 0 || n[2] >= dim[2] )
  {
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
    return;
  }
  
  double scale = min( min( voxel_size[0], voxel_size[1] ), voxel_size[2] ) * 1.0;

  vtkSmartPointer<vtkPolyDataAlgorithm> objsource;
  if ( GetProperties()->GetTensorRepresentation() == LayerPropertiesMRI::TR_Ellipsoid )
    objsource = vtkSmartPointer<vtkSphereSource>::New();
  else 
    objsource = vtkSmartPointer<vtkCubeSource>::New();
  objsource->Update();
  vtkPolyData* srcpolydata = objsource->GetOutput();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  scalars->SetNumberOfComponents( 4 );
//  srcpolydata->GetPointData()->SetNormals( NULL );    // remove normals
  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  double pt[3];
  int nSkip = 1;
  switch ( nPlane )
  {
    case 0:
      for ( int i = 0; i < dim[1]; i+=nSkip )
      {
        for ( int j = 0; j < dim[2]; j+=nSkip )
        {
          pt[0] = orig[0] + voxel_size[0] * n[0];
          pt[1] = orig[1] + voxel_size[1] * i;
          pt[2] = orig[2] + voxel_size[2] * j;           
          BuildTensorGlyph( imagedata, n[0], i, j, pt, scale, srcpolydata, scalars, append ) ;
        }
      }
      break;
    case 1:
      for ( int i = 0; i < dim[0]; i+=nSkip )
      {
        for ( int j = 0; j < dim[2]; j+=nSkip )
        {
          pt[0] = orig[0] + voxel_size[0] * i;
          pt[1] = orig[1] + voxel_size[1] * n[1];
          pt[2] = orig[2] + voxel_size[2] * j;            
          BuildTensorGlyph( imagedata, i, n[1], j, pt, scale, srcpolydata, scalars, append );
        }
      }
      break;
    case 2:
      for ( int i = 0; i < dim[0]; i+=nSkip )
      {
        for ( int j = 0; j < dim[1]; j+=nSkip )
        {
          pt[0] = orig[0] + voxel_size[0] * i;
          pt[1] = orig[1] + voxel_size[1] * j;
          pt[2] = orig[2] + voxel_size[2] * n[2];
          BuildTensorGlyph( imagedata, i, j, n[2], pt, scale, srcpolydata, scalars, append );
        }
      }
      break;
    default:
      break;
  }
  append->Update();
  vtkPolyData* polydata = append->GetOutput();
  polydata->GetPointData()->SetScalars( scalars );
  vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
  vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
}

void LayerMRI::BuildTensorGlyph( vtkImageData* imagedata,
                                         int i, int j, int k, 
                                         double* pt, double scale, 
                                         vtkPolyData* sourcepolydata,
                                         vtkUnsignedCharArray* scalars,
                                         vtkPolyDataAlgorithm* a)
{
  double** D = private_buf1_3x3;
  double** v = private_buf2_3x3;
  double w[3];
  D[0][0] = imagedata->GetScalarComponentAsDouble( i, j, k, 0 );
  D[0][1] = imagedata->GetScalarComponentAsDouble( i, j, k, 1 );
  D[0][2] = imagedata->GetScalarComponentAsDouble( i, j, k, 2 );
  D[1][0] = imagedata->GetScalarComponentAsDouble( i, j, k, 3 );
  D[1][1] = imagedata->GetScalarComponentAsDouble( i, j, k, 4 );
  D[1][2] = imagedata->GetScalarComponentAsDouble( i, j, k, 5 );
  D[2][0] = imagedata->GetScalarComponentAsDouble( i, j, k, 6 );
  D[2][1] = imagedata->GetScalarComponentAsDouble( i, j, k, 7 );
  D[2][2] = imagedata->GetScalarComponentAsDouble( i, j, k, 8 );
  if ( vtkMath::Jacobi( D, w, v ) )
  {
    v[1][0] = -v[1][0];         // by default invert Y !!
    v[1][1] = -v[1][1];
    v[1][2] = -v[1][2];
    if ( GetProperties()->GetTensorInversion() == LayerPropertiesMRI::VI_X )
    {
      v[0][0] = -v[0][0];
      v[0][1] = -v[0][1];
      v[0][2] = -v[0][2];
    }
    else if ( GetProperties()->GetTensorInversion() == LayerPropertiesMRI::VI_Y )
    {
      v[1][0] = -v[1][0];
      v[1][1] = -v[1][1];
      v[1][2] = -v[1][2];
    }
    else if ( GetProperties()->GetTensorInversion() == LayerPropertiesMRI::VI_Z )
    {
      v[2][0] = -v[2][0];
      v[2][1] = -v[2][1];
      v[2][2] = -v[2][2];
    }
    
    // make the vectors in right hand coordinate
    double v0[3] = { v[0][0], v[1][0], v[2][0] };
    double v1[3] = { v[0][1], v[1][1], v[2][1] };
    double v2[3];
    vtkMath::Cross( v0, v1, v2 );
    v[0][2] = v2[0];
    v[1][2] = v2[1];
    v[2][2] = v2[2];    
    
    vtkMath::Normalize( w );
 //   double w_sum = 1;//fabs(w[0]) + fabs(w[1]) + fabs(w[2]);
    w[0] = fabs(w[0]*scale);
    w[1] = fabs(w[1]*scale);
    w[2] = fabs(w[2]*scale);
    
    vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
    tr->Identity();
    tr->Translate( pt );
    double m[16];
    memset( m, 0, sizeof(double)*16 );
    m[15] = 1;
    
    m[0] = v[0][0];
    m[1] = v[0][1];
    m[2] = v[0][2];

    m[4] = v[1][0];
    m[5] = v[1][1];
    m[6] = v[1][2];

    
    m[8] = v[2][0];
    m[9] = v[2][1];
    m[10]= v[2][2];
 
    tr->Concatenate(m);
    tr->Scale( w );
//    tr->RotateZ( -90 );
    unsigned char c[4];
    c[0] = (int)(fabs( v[0][0] *255 ) );
    c[1] = (int)(fabs( v[1][0] *255 ) );
    c[2] = (int)(fabs( v[2][0] *255 ) );
    c[3] = 255;
    int nPts = sourcepolydata->GetPoints()->GetNumberOfPoints();
    for ( int i = 0; i < nPts; i++ )
      scalars->InsertNextTupleValue( c );
    
    vtkSmartPointer<vtkTransformPolyDataFilter> filter = 
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    filter->SetTransform( tr );
    filter->SetInput( sourcepolydata );

    a->AddInput( filter->GetOutput() );
  }
}

void LayerMRI::GetRASCenter( double* rasPt )
{
  MRI* mri = m_volumeSource->GetMRITarget();
  ::MRIvoxelToWorld( mri, 
                     mri->width / 2 - 0.5,
                     mri->height / 2 - 0.5,
                     mri->depth / 2 - 0.5,
                     &rasPt[0], &rasPt[1], &rasPt[2] ); 
}


void LayerMRI::UpdateVoxelValueRange( double dValue )
{
  mProperties->UpdateValueRange( dValue );
}

// Get voxel value range of a selected rectangle region defined by pt0, pt1
bool LayerMRI::GetVoxelValueRange( const double* pt0, const double* pt1, int nPlane, double* range_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  
  if ( nPlane < 0 || nPlane >= dim[nPlane] )
    return false;
  
  // find the index range of the selection
  int n0[3] = { 10000000, 10000000, 10000000 }, 
      n1[3] = { -10000000, -10000000, -10000000 };
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 ); 
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      p0 = max( 0, min( dim[i]-1, p0 ) );
      p1 = max( 0, min( dim[i]-1, p1 ) );
      n0[i] = min( p0, min( p1, n0[i] ) );
      n1[i] = max( p0, max( p1, n0[i] ) );

    }
  }
  
  int nActiveComp = GetActiveFrame();
  range_out[0] = m_imageData->GetScalarComponentAsDouble( n0[0], n0[1], n0[2], nActiveComp );
  range_out[1] = range_out[0];
  
  for ( int i = n0[0]; i <= n1[0]; i++ )
  {
    for ( int j = n0[1]; j <= n1[1]; j++ )
    {
      for ( int k = n0[2]; k <= n1[2]; k++ )
      {
        double value = m_imageData->GetScalarComponentAsDouble( i, j, k, nActiveComp );
        if ( range_out[0] > value )
          range_out[0] = value;
        if ( range_out[1] < value )
          range_out[1] = value;
      }
    }
  }
  
  return true;
}

// Get rectangle region stats
bool LayerMRI::GetVoxelStats( const double* pt0, const double* pt1, int nPlane, double* mean_out, double* sd_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  
  if ( nPlane < 0 || nPlane >= dim[nPlane] )
    return false;
  
  // find the index range of the selection
  int n0[3] = { 10000000, 10000000, 10000000 }, 
  n1[3] = { -10000000, -10000000, -10000000 };
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 ); 
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      p0 = max( 0, min( dim[i]-1, p0 ) );
      p1 = max( 0, min( dim[i]-1, p1 ) );
      n0[i] = min( p0, min( p1, n0[i] ) );
      n1[i] = max( p0, max( p1, n0[i] ) );
    }
  }
  
  int nActiveComp = GetActiveFrame();  
  double dMean = 0;
  int nCount = 0;
  for ( int i = n0[0]; i <= n1[0]; i++ )
  {
    for ( int j = n0[1]; j <= n1[1]; j++ )
    {
      for ( int k = n0[2]; k <= n1[2]; k++ )
      {
        dMean += m_imageData->GetScalarComponentAsDouble( i, j, k, nActiveComp );
        nCount++;
      }
    }
  }
  if ( nCount > 0 )
    *mean_out = dMean / nCount;
  
  if ( sd_out )
  {
    double sd = 0;
    for ( int i = n0[0]; i <= n1[0]; i++ )
    {
      for ( int j = n0[1]; j <= n1[1]; j++ )
      {
        for ( int k = n0[2]; k <= n1[2]; k++ )
        {
          double value = m_imageData->GetScalarComponentAsDouble( i, j, k, nActiveComp );
          sd += ( value-(*mean_out) ) * ( value-(*mean_out) );
        }
      }
    }
    if (nCount > 1)
      *sd_out = sqrt( sd / (nCount-1) );
    else
      *sd_out = 0;
  }
  
  return true;
}

void LayerMRI::ResetWindowLevel()
{
  double range[2];
  m_imageData->GetScalarRange( range );
  GetProperties()->SetMinMaxGrayscaleWindow( range[0], range[1] ); 
  GetProperties()->SetMinMaxJetScaleWindow( range[0], range[1] );
  GetProperties()->SetHeatScale( range[0], (range[0]+range[1])/2, range[1] );
}

int LayerMRI::GetDataType()
{
  return ( m_volumeSource ? m_volumeSource->GetDataType() : -1 );
}

COLOR_TABLE* LayerMRI::GetEmbeddedColorTable()
{
  return ( m_volumeSource ? m_volumeSource->GetEmbeddedColorTable(): NULL );
}

void LayerMRI::SnagToVoxelCenter( const double* pt_in, double* pt_out )
{
  if ( m_imageData == NULL )
  {
    pt_out[0] = pt_in[0];
    pt_out[1] = pt_in[1];
    pt_out[2] = pt_in[2];
    return;
  }

  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
    pt_out[i] = ( (int)( (pt_in[i] - orig[i])/vsize[i] + 0.5 ) ) * vsize[i] + orig[i];
}
