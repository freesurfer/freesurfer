/**
 * @file  LayerMRI.cpp
 * @brief Layer class for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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
#include "vtkMarchingSquares.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkSimpleLabelEdgeFilter.h"
#include "vtkImageResample.h"
#include "vtkPolyDataWriter.h"

#include "LayerPropertiesMRI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "MainWindow.h"
#include "RenderView.h"
#include "vtkMath.h"
#include "BuildContourThread.h"
#include "Contour2D.h"
#include "SurfaceRegion.h"
#include "SurfaceRegionGroups.h"

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define IMAGE_RESAMPLE_FACTOR     4.0     // must be multiples of 2

LayerMRI::LayerMRI( LayerMRI* ref ) : LayerVolumeBase(),
    m_volumeSource( NULL),
    m_volumeRef( ref ? ref->GetSourceVolume() : NULL ),
    m_bResampleToRAS( false ),
    m_bReorient( false ),
    m_nSampleMethod( SAMPLE_NEAREST ),
    m_bConform( false ),
    m_currentSurfaceRegion( NULL )
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

  m_actorContour = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
  m_actorContour->SetMapper( mapper );
  
  m_propVolume = vtkSmartPointer<vtkVolume>::New();
  
  m_nThreadID = 0;
  m_surfaceRegionGroups = new SurfaceRegionGroups( this );
  
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
	
  if ( m_sFilename.size() > 0 )
    GetProperties()->SaveSettings( m_sFilename.c_str() );
	
	if ( m_volumeSource )
		delete m_volumeSource;
    
  for ( size_t i = 0; i < m_segActors.size(); i++ )
  {
    m_segActors[i].actor->Delete();
  }
  
  delete m_surfaceRegionGroups;
  for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
  {
    delete m_surfaceRegions[i];
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

const char* LayerMRI::GetOrientationString()
{
  if ( m_volumeSource )
    return m_volumeSource->GetOrientationString();
  else
    return "RAS";
}

void LayerMRI::SetConform( bool bConform )
{
  m_bConform = bConform;
}

bool LayerMRI::LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_volumeSource )
    delete m_volumeSource;

  m_volumeSource = new FSVolume( m_volumeRef );
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_volumeSource->SetConform( m_bConform );
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

  GetProperties()->SetVolumeSource( m_volumeSource );
  GetProperties()->RestoreSettings( m_sFilename.c_str() );
  
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

    GetProperties()->SetVolumeSource( m_volumeSource );
    // mProperties->SetColorMap( LayerPropertiesMRI::LUT );

    m_sFilename = "";

    if ( bCopyVoxelData )
    {
      GetProperties()->CopySettings( mri->GetProperties() ); 
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
  
  if ( !m_volumeSource->UpdateMRIFromImage( m_imageData, wnd, event, !m_bReorient ) )
    return false;
  
  // now first save a copy of the old file if required
  if ( MainWindow::GetMainWindowPointer()->GetSaveCopy() )
  {
    wxString new_fn = wxString(m_sFilename.c_str()) + "~";
    rename( m_sFilename.c_str(), new_fn.c_str() );
  }
   
  bool bSaved = m_volumeSource->MRIWrite( m_sFilename.c_str(), 
                                          (mReslice[0]->GetInterpolationMode() == VTK_RESLICE_NEAREST ? SAMPLE_NEAREST : SAMPLE_TRILINEAR ) );
  m_bModified = !bSaved;

  event.SetInt( 99 );
  wxPostEvent( wnd, event );

  return bSaved;
}

bool LayerMRI::IsTransformed()
{
  vtkMatrix4x4* mat = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() )->GetMatrix();
  return !MyUtils::IsIdentity( mat->Element );
}

void LayerMRI::DoRestore()
{
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform(); 
  slice_tr->Identity();
  ras_tr->Identity();
  for ( int i = 0; i < 3; i++ )
    mReslice[i]->Modified();
}

bool LayerMRI::DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT || rotations[0].SampleMethod == SAMPLE_NEAREST )
      mReslice[i]->SetInterpolationModeToNearestNeighbor();
    else
      mReslice[i]->SetInterpolationModeToLinear();
  }
  
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  // also record transformation in RAS space
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  
  for ( size_t i = 0; i < rotations.size(); i++ )
  {
    double v[3] = { 0, 0, 0 };
    v[rotations[i].Plane] = 1;
    double dTargetPoint[3];
    RASToTarget( rotations[i].Point, dTargetPoint );
    
    slice_tr->Translate( dTargetPoint[0], dTargetPoint[1], dTargetPoint[2] );
    slice_tr->RotateWXYZ( -rotations[i].Angle, v );
    slice_tr->Translate( -dTargetPoint[0], -dTargetPoint[1], -dTargetPoint[2] );
    
    // record transformation in RAS space
    ras_tr->Translate( -rotations[i].Point[0], -rotations[i].Point[1], -rotations[i].Point[2] );
    double vp[3];
    vp[0] = dTargetPoint[0] + v[0];
    vp[1] = dTargetPoint[1] + v[1];
    vp[2] = dTargetPoint[2] + v[2];
    TargetToRAS( vp, vp );
    v[0] = vp[0] - rotations[i].Point[0];
    v[1] = vp[1] - rotations[i].Point[1];
    v[2] = vp[2] - rotations[i].Point[2];    
    ras_tr->RotateWXYZ( rotations[i].Angle, v );
    ras_tr->Translate( rotations[i].Point[0], rotations[i].Point[1], rotations[i].Point[2] );
  }

  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->Modified();
  }
  
  return true;
}

void LayerMRI::DoTranslate( double* offset )
{
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  slice_tr->Translate( -offset[0], -offset[1], -offset[2] );
  ras_tr->Translate( offset );
  for ( int i = 0; i < 3; i++ )
    mReslice[i]->Modified();
}


void LayerMRI::DoScale( double* scale, int nSampleMethod )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT || nSampleMethod == SAMPLE_NEAREST )
      mReslice[i]->SetInterpolationModeToNearestNeighbor();
    else
      mReslice[i]->SetInterpolationModeToLinear();
  }
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  double cpt[3], target_cpt[3];
  GetRASCenter( cpt );
  RASToTarget( cpt, target_cpt );
  slice_tr->Translate( target_cpt[0], target_cpt[1], target_cpt[2] );
  slice_tr->Scale( 1.0/scale[0], 1.0/scale[1], 1.0/scale[2] );
  slice_tr->Translate( -target_cpt[0], -target_cpt[1], -target_cpt[2] );
  
  ras_tr->Translate( -cpt[0], -cpt[1], -cpt[2] );
  ras_tr->Scale( scale );
  ras_tr->Translate( cpt[0], cpt[1], cpt[2] );
  
  for ( int i = 0; i < 3; i++ )
    mReslice[i]->Modified();
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
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->Identity();
  m_actorContour->SetUserTransform( tr->GetLinearInverse() );
  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
    mReslice[i]->SetInput( m_imageData );
    mReslice[i]->BorderOff();
    mReslice[i]->SetResliceTransform( tr );
//    mReslice[i]->AutoCropOutputOn();

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
    mColorMap[i]->SetLookupTable( GetProperties()->GetGrayScaleTable() );
    mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    mColorMap[i]->SetOutputFormatToRGBA();
    mColorMap[i]->PassAlphaToOutputOn();

    //
    // Prop in scene with plane mesh and texture.
    //
    m_sliceActor2D[i]->SetInput( mColorMap[i]->GetOutput() );
    if ( !MainWindow::GetMainWindowPointer()->GetRenderView(3)->GetRenderDisabled())
      m_sliceActor3D[i]->SetInput( mColorMap[i]->GetOutput() );

    mEdgeFilter[i] = vtkSmartPointer<vtkSimpleLabelEdgeFilter>::New();
    mResample[i] = vtkSmartPointer<vtkImageResample>::New();
    mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR );
    mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR );
    mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR );
    mResample[i]->SetInterpolationModeToNearestNeighbor();
    
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
    m_sliceActor2D[i]->SetOpacity( GetProperties()->GetOpacity() );
    m_sliceActor3D[i]->SetOpacity( GetProperties()->GetOpacity() );
  }
  m_actorContour->GetProperty()->SetOpacity( GetProperties()->GetOpacity() );
}

void LayerMRI::UpdateColorMap ()
{
  assert( GetProperties() );

  for ( int i = 0; i < 3; i++ )
    mColorMap[i]->SetActiveComponent( m_nActiveFrame );

  for ( int i = 0; i < 3; i++ )
    mColorMap[i]->SetLookupTable( GetProperties()->GetActiveLookupTable() );
  
  m_actorContour->GetMapper()->SetLookupTable( GetProperties()->GetActiveLookupTable() );
}

void LayerMRI::UpdateResliceInterpolation ()
{
  assert( GetProperties() );

  for ( int i = 0; i < 3; i++ )
  {
    if ( mReslice[i].GetPointer() )
    {
      mReslice[i]->SetInterpolationMode( GetProperties()->GetResliceInterpolation() );
    }
  }
}

void LayerMRI::UpdateTextureSmoothing ()
{
  assert( GetProperties() );

  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetInterpolate( GetProperties()->GetTextureSmoothing() );
    m_sliceActor3D[i]->SetInterpolate( GetProperties()->GetTextureSmoothing() );
  }
}

void LayerMRI::UpdateContour( int nSegValue )
{  
  UpdateContourActor( nSegValue );
}

void LayerMRI::UpdateContourActor( int nSegValue )
{
  // Generate a new thread id before creating the thread. so that mainwindow will be able to determine 
  // if a build contour result is already expired, by comparing the returned id and current id. If they
  // are different, it means a new thread is rebuilding the contour
  m_nThreadID++;
  BuildContourThread* thread = new BuildContourThread( MainWindow::GetMainWindowPointer() );
  thread->SetSmoothFactor( GetProperties()->GetContourSmoothIterations() );
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
  if ( GetProperties()->GetContourUseImageColorMap() )
  {
    m_actorContour->GetMapper()->ScalarVisibilityOn();
  }
  else
  {
    m_actorContour->GetMapper()->ScalarVisibilityOff();
    m_actorContour->GetProperty()->SetColor( GetProperties()->GetContourColor() );
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

  if ( GetProperties()->GetDisplayVector() || GetProperties()->GetDisplayTensor() )
    renderer->AddViewProp( m_glyphActor2D[nPlane] );
  else
    renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerMRI::Remove2DProps( vtkRenderer* renderer, int nPlane )
{
  wxASSERT ( nPlane >= 0 && nPlane <= 2 );

  if ( GetProperties()->GetDisplayVector() || GetProperties()->GetDisplayTensor() )
    renderer->RemoveViewProp( m_glyphActor2D[nPlane] );
  else
    renderer->RemoveViewProp( m_sliceActor2D[nPlane] );
}

void LayerMRI::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
	bool bContour = GetProperties()->GetShowAsContour();
	for ( int i = 0; i < 3; i++ )
	{
		if ( !bContour && ( bSliceVisibility == NULL || bSliceVisibility[i] ) )
    {
      if ( GetProperties()->GetDisplayVector() || GetProperties()->GetDisplayTensor() )
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
    for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
      m_surfaceRegions[i]->AppendProps( renderer );
	}
}

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

  assert( GetProperties() );
 
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
      m_sliceActor2D[0]->RotateY( 90 );
      m_sliceActor3D[0]->PokeMatrix( matrix );
      m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
      m_sliceActor3D[0]->RotateX( 90 );
      m_sliceActor3D[0]->RotateY( 90 );
  
      // Putting negatives in the reslice axes cosines will flip the
      // image on that axis.
      mReslice[0]->SetResliceAxesDirectionCosines( 0, 1, 0,
          0, 0, 1,
          1, 0, 0 );
      mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
      mReslice[0]->Modified();
      break;
    case 1:
      m_sliceActor2D[1]->PokeMatrix( matrix );
      m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
      m_sliceActor2D[1]->RotateX( 90 );
      m_sliceActor3D[1]->PokeMatrix( matrix );
      m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
      m_sliceActor3D[1]->RotateX( 90 );
  
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
  if ( GetProperties()->GetDisplayVector() )
  {
    UpdateVectorActor( nPlane );
  }
  else if ( GetProperties()->GetDisplayTensor() )
  {
    UpdateTensorActor( nPlane );
  }
  else if ( /*GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT &&*/ GetProperties()->GetShowLabelOutline() )
  {
    UpdateLabelOutline();
  }
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
      if ( GetProperties()->GetDisplayTensor() && 
           GetProperties()->GetTensorRepresentation() == LayerPropertiesMRI::TR_Ellipsoid )
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
    if ( GetProperties()->GetDisplayVector() )
      this->UpdateVectorActor();
    else if ( GetProperties()->GetDisplayTensor() )
      this->UpdateTensorActor();
    this->SendBroadcast( iMessage, this, this );
  }
  else if ( iMessage == "LayerLabelOutlineChanged" )
  {
    UpdateLabelOutline();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "LayerUpSampleMethodChanged" )
  {
    UpdateUpSampleMethod();
    this->SendBroadcast( "LayerActorUpdated", this, this );
  }
  else if ( iMessage == "SurfaceRegionColorChanged" )
    this->SendBroadcast( iMessage, this, this );
  
  LayerVolumeBase::DoListenToMessage( iMessage, iData, sender );
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

  vtkAbstractTransform* tr = mReslice[0]->GetResliceTransform();
  double pos_new[3];
  tr->TransformPoint( pos, pos_new );
  
  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  int* ext = m_imageData->GetExtent();
  
  
  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos_new[i] - orig[i] ) / vsize[i] + 0.5 );
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

void LayerMRI::NativeRASToTkReg( const double* pos_in, double* pos_out )
{
  m_volumeSource->NativeRASToTkReg( pos_in, pos_out );
}

void LayerMRI::TkRegToNativeRAS( const double* pos_in, double* pos_out )
{
  m_volumeSource->TkRegToNativeRAS( pos_in, pos_out );
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
  if ( GetProperties()->GetShowAsContour() )
  {
    if ( m_actorContour.GetPointer() == prop )
      return true;
    else
    {
      for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
      {
        if ( m_surfaceRegions[i]->GetMeshActor() == prop )
          return true;
      }
      return false;
    }
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      if ( m_sliceActor3D[i] == prop )
        return true;
    }
    return false;
  }
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
  GetProperties()->UpdateValueRange( dValue );
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
bool LayerMRI::GetVoxelStatsRectangle( const double* pt0, const double* pt1, int nPlane, double* mean_out, double* sd_out, int* cnt_out )
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
      n1[i] = max( p0, max( p1, n1[i] ) );
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
  
  if ( cnt_out )
    *cnt_out = nCount;
  
  return true;
}


// memory allocated for indice_out and value_out need to be freed outside of this function!
bool LayerMRI::GetVoxelsOnLine( const double* pt0, const double* pt1, int nPlane, int*& indice_out, double*& value_out, int* cnt_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  
  if ( nPlane < 0 || nPlane >= dim[nPlane] )
    return false;
  
  // find the index range of the selection
  int n0[3], n1[3]; 
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 ); 
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 ); 
      n0[i] = max( 0, min( dim[i]-1, p0 ) );
      n1[i] = max( 0, min( dim[i]-1, p1 ) );
    }
  }
  
  std::vector<int> indices = GetVoxelIndicesBetweenPoints( n0, n1 );
  std::vector<double> values;
  
  int nActiveComp = GetActiveFrame();  
  double dMean = 0;
  int nCount = 0;
  for ( size_t i = 0; i < indices.size(); i += 3 )
  {
    double value = m_imageData->GetScalarComponentAsDouble( indices[i], indices[i+1], indices[i+2], nActiveComp );
    dMean += value;
    values.push_back( value );
    nCount++;
  }
  
  indice_out = new int[nCount*3];
  value_out = new double[nCount];
  double pt[3], ras[3];
  int nIndex[3];
  for ( int i = 0; i < nCount; i++ )
  {
    pt[0] = indices[i*3]*voxel_size[0] + orig[0];
    pt[1] = indices[i*3+1]*voxel_size[1] + orig[1];
    pt[2] = indices[i*3+2]*voxel_size[2] + orig[2];
    RemapPositionToRealRAS( pt, ras );
    RASToOriginalIndex( ras, nIndex );
    indice_out[i*3]   = nIndex[0];
    indice_out[i*3+1] = nIndex[1];
    indice_out[i*3+2] = nIndex[2];
    value_out[i] = values[i];
  }
  
  *cnt_out = nCount;
  return true;
}

std::vector<int> LayerMRI::GetVoxelIndicesBetweenPoints( int* n0, int* n1 )
{
  std::vector<int> indices;
  if ( n1[0] == n0[0] && n1[1] == n0[1] && n1[2] == n0[2] )
  {
    indices.push_back( n0[0] );
    indices.push_back( n0[1] );
    indices.push_back( n0[2] );
  }
  else if ( fabs(n0[0]-n1[0]) <= 1 && fabs(n0[1]-n1[1]) <= 1 && fabs(n0[2]-n1[2]) <= 1 )
  {
    indices.push_back( n0[0] );
    indices.push_back( n0[1] );
    indices.push_back( n0[2] );
    indices.push_back( n1[0] );
    indices.push_back( n1[1] );
    indices.push_back( n1[2] );
  }
  else
  {
    int n[3];
    for ( int i = 0; i < 3; i++ )
      n[i] = (int)( (n0[i]+n1[i]) / 2.0 + 0.5 );
    
    indices = GetVoxelIndicesBetweenPoints( n0, n );
    std::vector<int> indices1 = GetVoxelIndicesBetweenPoints( n, n1 );
    for ( size_t i = 3; i < indices1.size(); i++ )
      indices.push_back( indices1[i] );
  }
    
  return indices;
}

void LayerMRI::ResetWindowLevel()
{
  double range[2];
  m_imageData->GetScalarRange( range );
  GetProperties()->SetMinMaxGrayscaleWindow( range[0], range[1] ); 
  GetProperties()->SetMinMaxGenericThreshold( range[0], range[1] );
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

void LayerMRI::SnapToVoxelCenter( const double* pt_in, double* pt_out )
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

void LayerMRI::UpdateLabelOutline()
{
  if ( /*GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT &&*/ GetProperties()->GetShowLabelOutline() )
  {
    double* vsize = m_imageData->GetSpacing();
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR );
      mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR );
      mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR );
      mResample[i]->SetInterpolationModeToNearestNeighbor();
      double pos[3] = { vsize[0]/IMAGE_RESAMPLE_FACTOR/2, vsize[1]/IMAGE_RESAMPLE_FACTOR/2, vsize[2]/IMAGE_RESAMPLE_FACTOR/2 };
      mResample[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      mEdgeFilter[i]->SetInputConnection( mResample[i]->GetOutputPort() );
      mColorMap[i]->SetInputConnection( mEdgeFilter[i]->GetOutputPort() );
      pos[i] = m_dSlicePosition[i];
      m_sliceActor2D[i]->SetPosition( pos );
      m_sliceActor3D[i]->SetPosition( pos );
    }
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      double pos[3] = { 0, 0, 0};
      mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      pos[i] = m_dSlicePosition[i]; 
      m_sliceActor2D[i]->SetPosition( pos );
      m_sliceActor3D[i]->SetPosition( pos );
    }
  }
  this->SendBroadcast( "ResampleFactorChanged", this );
}

void LayerMRI::UpdateUpSampleMethod()
{
  switch ( GetProperties()->GetUpSampleMethod() )
  {
    case LayerPropertiesMRI::UM_NearestNeighbor:
      for ( int i = 0; i < 3; i++ )
        mResample[i]->SetInterpolationModeToNearestNeighbor();
      break;
    case LayerPropertiesMRI::UM_BiLinear:
      for ( int i = 0; i < 3; i++ )
        mResample[i]->SetInterpolationModeToLinear();
      break;
    default:
      break;
  }
  if ( GetProperties()->GetUpSampleMethod() == LayerPropertiesMRI::UM_None )
  {
    for ( int i = 0; i < 3; i++ )
    {
      mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    }
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      mColorMap[i]->SetInputConnection( mResample[i]->GetOutputPort() );
      if ( !GetProperties()->GetShowLabelOutline() )
      {
        mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR/2 );
        mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR/2 );
        mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR/2 );
      }
    }
  }
  
  this->SendBroadcast( "ResampleFactorChanged", this );
}


void LayerMRI::GetCurrentLabelStats( int nPlane, float* label_out, int* count_out, float* area_out )
{
  if ( !m_imageData || nPlane < 0 || nPlane > 2 )
    return;
  
  double* origin = m_imageData->GetOrigin();
  int* dim = m_imageData->GetDimensions();
  double* pos = GetSlicePosition();
  double vs[3];
  m_imageData->GetSpacing( vs );
  
  int n[3];
  for ( int i = 0; i < 3; i++ )
    n[i] = (int)( ( pos[i] - origin[i] ) / vs[i] );
  
  float fLabel = 0;
  if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] && n[2] >= 0 && n[2] < dim[2] )
    fLabel = m_imageData->GetScalarComponentAsFloat( n[0], n[1], n[2], m_nActiveFrame );
  
  int cnt = 0;
  int ext[3][2] = { { 0, dim[0]-1 }, {0, dim[1]-1}, {0, dim[2]-1} };
  ext[nPlane][0] = ext[nPlane][1] = n[nPlane];
  for ( int i = ext[0][0]; i <= ext[0][1]; i++ )
  {
    for ( int j = ext[1][0]; j <= ext[1][1]; j++ )
    {
      for ( int k = ext[2][0]; k <= ext[2][1]; k++ )
      {
        if ( m_imageData->GetScalarComponentAsFloat( i, j, k, m_nActiveFrame ) == fLabel )
        {
          cnt++;
        }
      }
    }
  }
  vs[nPlane] = 1.0;
  
  *label_out = fLabel;
  *count_out = cnt;
  *area_out = cnt*vs[0]*vs[1]*vs[2]; 
}

vtkImageData* LayerMRI::GetSliceImageData( int nPlane )
{
  return mReslice[nPlane]->GetOutput();
}

bool LayerMRI::FloodFillByContour2D( double* ras, Contour2D* c2d )
{
  int nPlane = c2d->GetPlane();
  vtkImageData* image = c2d->GetThresholdedImage();
  vtkImageData* original_image = GetSliceImageData( nPlane );
  int* nDim = image->GetDimensions();      // 2D image
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  int nx = nDim[0], ny = nDim[1], x = 0, y = 0;
  switch ( nPlane )
  {
    case 0:
      x = n[1];
      y = n[2];
      break;
    case 1:
      x = n[0];
      y = n[2];
      break;
    case 2:
      x = n[0];
      y = n[1];
      break;
  }
  
  char* mask = new char[nDim[0]*nDim[1]];
  memset( mask, 0, nDim[0]*nDim[1] );
  double dVoxelValue = original_image->GetScalarComponentAsDouble( x, y, 0, 0 );
  double dMaskValue = image->GetScalarComponentAsDouble( x, y, 0, 0 );
  for ( int i = 0; i < nDim[0]; i++ )
  {
    for ( int j = 0; j < nDim[1]; j++ )
    {
      if ( original_image->GetScalarComponentAsDouble( i, j, 0, 0 ) == dVoxelValue &&
           image->GetScalarComponentAsDouble( i, j, 0, 0 ) == dMaskValue )
        mask[j*nDim[0]+i] = 1;
    }
  }
  
  int nFillValue = 2;  
  MyUtils::FloodFill( mask, x, y, nx, ny, nFillValue, 0 ); 
  
  int nActiveComp = this->GetActiveFrame();
  int cnt = 0;
  switch ( nPlane )
  {
    case 0:
      for ( int i = 0; i < nx; i++ )
      {
        for ( int j = 0; j < ny; j++ )
        {
          if ( mask[j*nx+i] == nFillValue )
          {
            m_imageData->SetScalarComponentFromFloat( n[nPlane], i, j, nActiveComp, m_fFillValue );
            cnt++;
          }
        }
      }
      break;
    case 1:
      for ( int i = 0; i < nx; i++ )
      {
        for ( int j = 0; j < ny; j++ )
        {
          if ( mask[j*nx+i] == nFillValue )
          {
            m_imageData->SetScalarComponentFromFloat( i, n[nPlane], j, nActiveComp, m_fFillValue );
            cnt++;
          }
        }
      }
      break;
    case 2:
      for ( int i = 0; i < nx; i++ )
      {
        for ( int j = 0; j < ny; j++ )
        {
          if ( mask[j*nx+i] == nFillValue )
          {
            m_imageData->SetScalarComponentFromFloat( i, j, n[nPlane], nActiveComp, m_fFillValue );
            cnt++;
          }
        }
      }
      break;
  }
  SetModified();
  this->SendBroadcast( "LayerActorUpdated", this );
  return true;
}

SurfaceRegion* LayerMRI::CreateNewSurfaceRegion( double* pt )
{
  SurfaceRegion* r = new SurfaceRegion( this );
  r->SetInput( vtkPolyData::SafeDownCast( m_actorContour->GetMapper()->GetInput() ) );
  r->AddPoint( pt );
  r->SetId( m_surfaceRegions.size() + 1 );
  m_surfaceRegions.push_back( r );
  int nGroup = 1;
  if ( m_currentSurfaceRegion )
  {
    m_currentSurfaceRegion->Highlight( false );
    nGroup = m_currentSurfaceRegion->GetGroup();
  }
  m_currentSurfaceRegion = r;
  r->SetGroup( nGroup );
  this->SendBroadcast( "SurfaceRegionAdded", this );
  r->AddListener( this );
  return r;
}
  
void LayerMRI::AddSurfaceRegionLoopPoint( double* pt )
{
  if ( m_currentSurfaceRegion )
  {
    m_currentSurfaceRegion->AddPoint( pt );
    this->SendBroadcast( "SurfaceRegionUpdated", this );
  }
}

void LayerMRI::CloseSurfaceRegion()
{
  if ( m_currentSurfaceRegion )
  {
    if ( !m_currentSurfaceRegion->Close() )
      DeleteCurrentSurfaceRegion();
    this->SendBroadcast( "SurfaceRegionUpdated", this );    
  }
}

SurfaceRegion* LayerMRI::SelectSurfaceRegion( double* pos )
{
  for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( m_surfaceRegions[i]->HasPoint( pos ) )
    {
      if ( m_currentSurfaceRegion )
        m_currentSurfaceRegion->Highlight( false );
      m_currentSurfaceRegion = m_surfaceRegions[i];
      m_currentSurfaceRegion->Highlight( true );
      return m_currentSurfaceRegion;
    }
  }
  
  return NULL;
}

SurfaceRegion* LayerMRI::SelectSurfaceRegion( int nId )
{
  for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( m_surfaceRegions[i]->GetId() == nId )
    {
      if ( m_currentSurfaceRegion )
        m_currentSurfaceRegion->Highlight( false );
      m_currentSurfaceRegion = m_surfaceRegions[i];
      m_currentSurfaceRegion->Highlight( true );
      return m_currentSurfaceRegion;
    }
  }
  return NULL;
}

bool LayerMRI::DeleteCurrentSurfaceRegion()
{
  if ( m_currentSurfaceRegion )
  {
    for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
    {
      if ( m_surfaceRegions[i] == m_currentSurfaceRegion )
      {
        m_surfaceRegions.erase( m_surfaceRegions.begin() + i );
        delete m_currentSurfaceRegion;
        m_currentSurfaceRegion = NULL;
        ResetSurfaceRegionIds();
        this->SendBroadcast( "SurfaceRegionRemoved", this );  
        return true;
      }
    }
  }
  return false;
}

void LayerMRI::ResetSurfaceRegionIds()
{
  for ( int i = 0; i < (int)m_surfaceRegions.size(); i++ )
  {
    m_surfaceRegions[i]->SetId( i+1 );
  }
}

bool LayerMRI::SaveAllSurfaceRegions( wxString& fn )
{
  FILE* fp = fopen( fn.c_str(), "w" );
  if ( !fp )
    return false;
  
  bool ret = SurfaceRegion::WriteHeader( fp, this, m_surfaceRegions.size() );
  for ( size_t i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( !m_surfaceRegions[i]->WriteBody( fp ) )
      ret = false;
  } 
  fclose( fp ); 
  return ret;
}

bool LayerMRI::LoadSurfaceRegions( wxString& fn )
{
  FILE* fp = fopen( fn.c_str(), "r" );
  if ( !fp )
  {
    cerr << "Can not open file " << fn.c_str() << endl;
    return false;
  }
  int nNum = 0;
  char ch[1000];
  float dTh_low, dTh_high;
  fscanf( fp, "VOLUME_PATH %s\nVOLUME_THRESHOLD %f %f\nSURFACE_REGIONS %d", ch, &dTh_low, &dTh_high, &nNum );
  if ( nNum == 0 )
    return false;
  
  bool bSuccess = true;
  for ( int i = 0; i < nNum; i++ )
  {
    SurfaceRegion* r = new SurfaceRegion( this );
    if ( r->Load( fp ) )
    {
      r->SetInput( vtkPolyData::SafeDownCast( m_actorContour->GetMapper()->GetInput() ) );
      m_surfaceRegions.push_back( r );
      r->Highlight( false );
      r->AddListener( this );
    }
    else
    {
      fclose( fp );
      bSuccess = false;
      break;
    }
  }
  
  ResetSurfaceRegionIds();
  if ( bSuccess )
  {
    GetProperties()->SetContourThreshold( dTh_low, dTh_high );
    GetProperties()->SetShowAsContour( true );
  }
  
  this->SendBroadcast( "SurfaceRegionAdded", this );
  return bSuccess;
}

void LayerMRI::SetCroppingBounds( double* bounds )
{
  m_volumeSource->SetCroppingBounds( bounds );
}

void LayerMRI::GetDisplayBounds( double* bounds )
{
  m_imageData->GetBounds( bounds );
  if ( mReslice[0].GetPointer() && mReslice[0]->GetAutoCropOutput() )
  {
    double d[6];
    m_imageData->GetBounds( d );
    vtkTransform* tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
    double pt[3];
    for ( int i = 0; i < 2; i++ )
    {
      for ( int j = 0; j < 2; j++ )
      {
        for ( int k = 0; k < 2; k++ )
        {
          pt[0] = d[i];
          pt[1] = d[2+j];
          pt[2] = d[4+k];
          tr->GetLinearInverse()->TransformPoint( pt, pt );
          if ( pt[0] < bounds[0] ) bounds[0] = pt[0];
          if ( pt[0] > bounds[1] ) bounds[1] = pt[0];
          if ( pt[1] < bounds[2] ) bounds[2] = pt[1];
          if ( pt[1] > bounds[3] ) bounds[3] = pt[1];
          if ( pt[2] < bounds[4] ) bounds[4] = pt[2];
          if ( pt[2] > bounds[5] ) bounds[5] = pt[2];
        }
      }
    }
  }
}

bool LayerMRI::SaveRegistration( const char* filename )
{
  return m_volumeSource->SaveRegistration( filename );
}
  
struct LabelStatsPrivate
{
  int id;
  int count;
  std::vector<double> values;
};

void LayerMRI::GetLabelStats( LayerMRI* label, int nPlane,
                              std::vector<int>& ids, 
                              std::vector<int>& numbers, 
                              std::vector<double>& means, 
                              std::vector<double>& stds )
{
  vtkImageData* mri_image = mReslice[nPlane]->GetOutput();
  vtkImageData* label_image = label->mReslice[nPlane]->GetOutput();
  int*    mri_dim = mri_image->GetDimensions();
  double* mri_vs = mri_image->GetSpacing();
  double* mri_orig = mri_image->GetOrigin();
  int*    label_dim = label_image->GetDimensions();
  double* label_vs = label_image->GetSpacing();
  double* label_orig = label_image->GetOrigin();
  
  // first find all label ids
  std::vector<LabelStatsPrivate> labels;
  for ( int i = 0; i < label_dim[0]; i++ )
  {
    for ( int j = 0; j < label_dim[1]; j++ )
    {
      int nId = (int)label_image->GetScalarComponentAsDouble( i, j, 0, 0 );
      if ( nId > 0 )
      {
        int mi = (int)( ( i*label_vs[0] + label_orig[0] - mri_orig[0] ) / mri_vs[0] );
        int mj = (int)( ( j*label_vs[1] + label_orig[1] - mri_orig[1] ) / mri_vs[1] );
        if ( mi >= 0 && mi < mri_dim[0] && mj >= 0 && mj < mri_dim[1] )
        {
          bool bFound = false;
          for ( size_t n = 0; n < labels.size(); n++ )
          {
            if ( nId == labels[n].id )
            {
              labels[n].values.push_back( mri_image->GetScalarComponentAsDouble( mi, mj, 0, 0 ) );
              labels[n].count++;
              bFound = true;
              break;
            }
            else if ( nId < labels[n].id )
            {
              LabelStatsPrivate l;
              l.id = nId;
              l.count = 1;
              l.values.push_back( mri_image->GetScalarComponentAsDouble( mi, mj, 0, 0 ) );
              labels.insert( labels.begin() + n, l );
              bFound = true;
              break; 
            }
          }
          if ( !bFound )
          {
            LabelStatsPrivate l;
            l.id = nId;
            l.count = 1;
            l.values.push_back( mri_image->GetScalarComponentAsDouble( mi, mj, 0, 0 ) );
            labels.push_back( l );
          }
        }
      }
    }
  }
  for ( size_t n = 0; n < labels.size(); n++ )
  {
    if ( labels[n].id > 0 )
    {
      ids.push_back( labels[n].id );
      numbers.push_back( labels[n].count );
      double dvalue = 0;
      for ( int i = 0; i < labels[n].count; i++ )
      {
        dvalue += labels[n].values[i];
      }
      double dmean = dvalue / labels[n].count;
      means.push_back( dmean );
      double sd = 0;
      for ( int i = 0; i < labels[n].count; i++ )
      {
        sd += (labels[n].values[i]-dmean)*(labels[n].values[i]-dmean);
      }
      if ( labels[n].count > 1 )
        sd = sqrt( sd / ( labels[n].count-1 ) );
      else
        sd = 0;  
      stds.push_back( sd );     
    }
  }
}

bool LayerMRI::SaveContourToFile( const char* filename )
{
  MATRIX* mat = m_volumeSource->GetTargetToRASMatrix();
  double m[16];
  for ( int i = 0; i < 16; i++ )
  {
    m[i] = (double) *MATRIX_RELT((mat),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &mat );
  vtkSmartPointer<vtkMatrix4x4> vmat = vtkSmartPointer<vtkMatrix4x4>::New();
  vmat->DeepCopy( m );
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->SetMatrix( vmat );
  vtkSmartPointer<vtkTransformPolyDataFilter> filter = 
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  filter->SetTransform( tr );
  filter->SetInput( vtkPolyDataMapper::SafeDownCast( m_actorContour->GetMapper())->GetInput() );
  filter->Update();
  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( filename );
  bool ret = writer->Write();
  writer->Delete();
  return ret;
}
