/**
 * @file  LayerSurface.cpp
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/03/20 19:03:53 $
 *    $Revision: 1.20 $
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
#include "LayerSurface.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkLODActor.h"
#include "vtkRGBATransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkDecimatePro.h"
#include "vtkMapperCollection.h"
#include "vtkTubeFilter.h"
#include "vtkPointData.h"
#include "LayerPropertiesSurface.h"
#include "MyUtils.h"
#include "FSSurface.h"
#include "LayerMRI.h"
#include "SurfaceOverlay.h"
#include "SurfaceOverlayProperties.h"

LayerSurface::LayerSurface( LayerMRI* ref ) : Layer(),
    m_surfaceSource( NULL ),
    m_bResampleToRAS( true ),
    m_volumeRef( ref ),
    m_nActiveOverlay( -1 )
{
  m_strTypeNames.push_back( "Surface" );

  for ( int i = 0; i < 3; i++ )
  {
    // m_nSliceNumber[i] = 0;
    m_sliceActor2D[i] = vtkActor::New();
    m_sliceActor3D[i] = vtkActor::New();
  }
  mProperties = new LayerPropertiesSurface();
  mProperties->AddListener( this );

// m_mainActor = vtkLODActor::New();
  m_mainActor = vtkActor::New();
  mLowResFilter = vtkSmartPointer<vtkDecimatePro>::New();
  mLowResFilter->SetTargetReduction( 0.9 );
// mMediumResFilter = vtkSmartPointer<vtkDecimatePro>::New();
// mMediumResFilter->SetTargetReduction( 0.9 );

  m_vectorActor = vtkActor::New();
  m_vectorActor->GetProperty()->SetColor( mProperties->GetVectorColor() );
  m_vectorActor->GetProperty()->SetPointSize( mProperties->GetVectorPointSize() );
}

LayerSurface::~LayerSurface()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->Delete();
    m_sliceActor3D[i]->Delete();
  }
  m_mainActor->Delete();
  m_vectorActor->Delete();

  delete mProperties;

  if ( m_surfaceSource )
    delete m_surfaceSource;
  
  for ( size_t i = 0; i < m_overlays.size(); i++ )
  {
    delete m_overlays[i];
  }
}

bool LayerSurface::LoadSurfaceFromFile( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_surfaceSource )
    delete m_surfaceSource;

  m_surfaceSource = new FSSurface( m_volumeRef ? m_volumeRef->GetSourceVolume() : NULL );
// m_surfaceSource->SetResampleToRAS( m_bResampleToRAS );
  if ( !m_surfaceSource->MRISRead( m_sFilename.c_str(), wnd, event,
                                   m_sVectorFilename.size() > 0 ? m_sVectorFilename.c_str() : NULL ) )
    return false;

  InitializeSurface();
  InitializeActors();

  mProperties->SetSurfaceSource( m_surfaceSource );

  event.SetInt( 100 );
  wxPostEvent( wnd, event );

  return true;
}


bool LayerSurface::LoadVectorFromFile( wxWindow* wnd, wxCommandEvent& event )
{
  if ( m_sVectorFilename.size() == 0 || !m_surfaceSource->MRISReadVectors( m_sVectorFilename.c_str(), wnd, event ) )
    return false;

  event.SetInt( 100 );
  wxPostEvent( wnd, event );

  this->SendBroadcast( "LayerModified", this );
  this->SendBroadcast( "LayerActorUpdated", this );

  return true;
}

bool LayerSurface::LoadCurvatureFromFile( const char* filename )
{
  if ( !m_surfaceSource->LoadCurvature( filename ) )
    return false;
  
  this->SendBroadcast( "LayerModified", this );
  this->SendBroadcast( "LayerCurvatureLoaded", this );
  return true;
}

bool LayerSurface::LoadOverlayFromFile( const char* filename )
{
  if ( !m_surfaceSource->LoadOverlay( filename ) )
    return false;
  
  // create overlay
  SurfaceOverlay* overlay = new SurfaceOverlay( this ); 
  std::string fn = filename;
  overlay->SetName( fn.substr( fn.find_last_of("/\\")+1 ).c_str() );
  m_overlays.push_back( overlay );
  
  SetActiveOverlay( m_overlays.size() - 1 );
  
  this->SendBroadcast( "LayerModified", this );
  this->SendBroadcast( "LayerOverlayAdded", this );
  return true; 
}

void LayerSurface::InitializeSurface()
{
  if ( m_surfaceSource == NULL )
    return;

  FSSurface* source = m_surfaceSource;

  float RASBounds[6];
  source->GetBounds( RASBounds );
  m_dWorldOrigin[0] = RASBounds[0];
  m_dWorldOrigin[1] = RASBounds[2];
  m_dWorldOrigin[2] = RASBounds[4];

  m_dWorldSize[0] = RASBounds[1] - RASBounds[0];
  m_dWorldSize[1] = RASBounds[3] - RASBounds[2];
  m_dWorldSize[2] = RASBounds[5] - RASBounds[4];
}


void LayerSurface::InitializeActors()
{
  if ( m_surfaceSource == NULL )
    return;

  // main surface actor
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( m_surfaceSource->GetPolyData() );
  m_mainActor->SetMapper( mapper );
  mapper->Update();

  // vector actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  /* vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
   tube->SetInput( m_surfaceSource->GetVectorPolyData() );
   tube->SetNumberOfSides( 6 );
   tube->SetRadius(0.05 );*/
  mapper->SetInput(  m_surfaceSource->GetVectorPolyData() );
  m_vectorActor->SetMapper( mapper );
  mapper->Update();
  /* mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
   mLowResFilter->SetInput( m_surfaceSource->GetPolyData() );
   mapper->SetInput( mLowResFilter->GetOutput() );
   m_mainActor->AddLODMapper( mapper );
   mapper->Update();
   mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
   mMediumResFilter->SetInput( m_surfaceSource->GetPolyData() );
   mapper->SetInput( mMediumResFilter->GetOutput() );
   m_mainActor->AddLODMapper( mapper );
   mapper->Update();*/

  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslicePlane[i] = vtkSmartPointer<vtkPlane>::New();
    mReslicePlane[i]->SetOrigin( 0, 0, 0 );
    mReslicePlane[i]->SetNormal( (i==0), (i==1), (i==2) );

    vtkSmartPointer<vtkCutter> cutter =
      vtkSmartPointer<vtkCutter>::New();
    cutter->SetInput( m_surfaceSource->GetPolyData() );
    cutter->SetCutFunction( mReslicePlane[i] );

    //
    // Mappers for the lines.
    //
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput( cutter->GetOutput() );
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInput( cutter->GetOutput() );
    //
    // Actors in the scene, drawing the mapped lines.
    //
    m_sliceActor2D[i]->SetMapper( mapper );
//  m_sliceActor2D[i]->SetBackfaceProperty( m_sliceActor2D[i]->MakeProperty() );
//  m_sliceActor2D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor2D[i]->SetProperty( m_sliceActor2D[i]->MakeProperty() );
    m_sliceActor2D[i]->GetProperty()->SetInterpolationToFlat();
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( mProperties->GetEdgeThickness() );

    m_sliceActor3D[i]->SetMapper( mapper2 );
//  m_sliceActor3D[i]->SetBackfaceProperty( m_sliceActor3D[i]->MakeProperty() );
//  m_sliceActor3D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor3D[i]->SetProperty( m_sliceActor3D[i]->MakeProperty() );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( mProperties->GetEdgeThickness() );
    m_sliceActor2D[i]->GetProperty()->SetInterpolationToFlat();

    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }

  this->UpdateOpacity();
  this->UpdateColorMap();
}

void LayerSurface::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    // m_sliceActor2D[i]->GetProperty()->SetOpacity( mProperties->GetOpacity() );
    // m_sliceActor3D[i]->SetOpacity( mProperties->GetOpacity() );
  }
  m_mainActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
}

void LayerSurface::UpdateEdgeThickness()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( mProperties->GetEdgeThickness() );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( mProperties->GetEdgeThickness() );
  }
  
  if ( mProperties->GetEdgeThickness() == 0 )
  {
    if ( IsVisible() )
    {
      for ( int i = 0; i < 3; i++ )
      {
        m_sliceActor2D[i]->SetVisibility( 0 );
        m_sliceActor3D[i]->SetVisibility( 0 );
      }
    }
  }
  else if ( IsVisible() && !m_sliceActor2D[0]->GetVisibility() )
  {
    for ( int i = 0; i < 3; i++ )
    {
      m_sliceActor2D[i]->SetVisibility( 1 );
      m_sliceActor3D[i]->SetVisibility( 1 );
    }
  }
}

void LayerSurface::UpdateVectorPointSize()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor->GetProperty()->SetPointSize( mProperties->GetVectorPointSize() );
  }
}

void LayerSurface::UpdateColorMap()
{
  if ( m_surfaceSource == NULL )
    return;

  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetProperty()->SetColor( mProperties->GetEdgeColor() );
    m_sliceActor3D[i]->GetProperty()->SetColor( mProperties->GetEdgeColor() );
    m_sliceActor2D[i]->GetMapper()->ScalarVisibilityOff();
    m_sliceActor3D[i]->GetMapper()->ScalarVisibilityOff();
  }

  m_mainActor->GetProperty()->SetColor( mProperties->GetBinaryColor() );
  m_vectorActor->GetProperty()->SetColor( mProperties->GetVectorColor() );
  if ( m_surfaceSource->IsCurvatureLoaded() )
  {
    if ( mProperties->GetCurvatureLUT() != m_mainActor->GetMapper()->GetLookupTable() )
      m_mainActor->GetMapper()->SetLookupTable( mProperties->GetCurvatureLUT() );
    /*  vtkSmartPointer<vtkMapperCollection> mc = m_mainActor->GetLODMappers();
      mc->InitTraversal();
      vtkMapper* mapper = NULL;
      while ( ( mapper = mc->GetNextItem() ) != NULL )
      {
       mapper->SetLookupTable( mProperties->GetCurvatureLUT() );
      } */
  }
  
  UpdateOverlay();
}

void LayerSurface::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  wxASSERT ( nPlane >= 0 && nPlane <= 2 );

  renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerSurface::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  renderer->AddViewProp( m_mainActor );
  renderer->AddViewProp( m_vectorActor );

  for ( int i = 0; i < 3; i++ )
    renderer->AddViewProp( m_sliceActor3D[i] );
}

/*
void LayerSurface::SetSliceNumber( int* sliceNumber )
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

void LayerSurface::SetSlicePositionToWorldCenter()
{
  if ( m_surfaceSource == NULL )
    return;

  // Get some values from the MRI.
  double pos[3];
  for ( int i = 0; i < 3; i++ )
    pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];

  SetSlicePosition( pos );
}

void LayerSurface::OnSlicePositionChanged( int nPlane )
{
  if ( m_surfaceSource == NULL )
    return;

  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  switch ( nPlane )
  {
  case 0:
    mReslicePlane[0]->SetOrigin( m_dSlicePosition[0], 0, 0  );
    m_sliceActor2D[0]->SetPosition( 0.1, 0, 0 );
    break;
  case 1:
    mReslicePlane[1]->SetOrigin( 0, m_dSlicePosition[1], 0 );
    m_sliceActor2D[1]->SetPosition( 0, 0.1, 0 );
    break;
  case 2:
    mReslicePlane[2]->SetOrigin( 0, 0, m_dSlicePosition[2]  );
    m_sliceActor2D[2]->SetPosition( 0, 0, -0.1 );
    break;
  }
}

LayerPropertiesSurface* LayerSurface::GetProperties()
{
  return mProperties;
}

void LayerSurface::DoListenToMessage( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->UpdateColorMap();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "OpacityChanged" )
  {
    this->UpdateOpacity();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "EdgeThicknessChanged" )
  {
    this->UpdateEdgeThickness();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "VectorPointSizeChanged" )
  {
    this->UpdateVectorPointSize();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
}

void LayerSurface::SetVisible( bool bVisible )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetVisibility( ( bVisible && mProperties->GetEdgeThickness() > 0 ) ? 1 : 0 );
    m_sliceActor3D[i]->SetVisibility( ( bVisible && mProperties->GetEdgeThickness() > 0 ) ? 1 : 0 );
  }
  m_mainActor->SetVisibility( bVisible ? 1 : 0 );
  m_vectorActor->SetVisibility( ( bVisible && m_surfaceSource && m_surfaceSource->GetActiveVector() >= 0 )? 1 : 0 );
  this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerSurface::IsVisible()
{
  return m_mainActor->GetVisibility() > 0;
}

bool LayerSurface::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( m_sliceActor2D[i] == prop || m_sliceActor3D[i] == prop )
      return true;
  }
  return prop == m_mainActor;
}

int LayerSurface::GetVertexIndexAtRAS( double* ras, double* distance )
{
  if ( m_surfaceSource == NULL )
    return -1;

  return m_surfaceSource->FindVertexAtRAS( ras, distance );
}

int LayerSurface::GetVertexIndexAtTarget( double* ras, double* distance )
{
  if ( m_surfaceSource == NULL )
    return -1;

  if ( m_volumeRef )
  {
    double realRas[3];
    m_volumeRef->TargetToRAS( ras, realRas );
    return m_surfaceSource->FindVertexAtRAS( realRas, distance );
  }
  else
    return m_surfaceSource->FindVertexAtRAS( ras, distance );
}


bool LayerSurface::GetRASAtVertex( int nVertex, double* ras )
{
  if ( m_surfaceSource == NULL )
    return false;

  return m_surfaceSource->GetRASAtVertex( nVertex, ras );
}


bool LayerSurface::GetTargetAtVertex( int nVertex, double* ras )
{
  if ( m_surfaceSource == NULL )
    return false;

  bool bRet = m_surfaceSource->GetRASAtVertex( nVertex, ras );
  if ( bRet && m_volumeRef )
    m_volumeRef->RASToTarget( ras, ras );

  return bRet;
}

void LayerSurface::SetActiveSurface( int nSurface )
{
  if ( m_surfaceSource && m_surfaceSource->SetActiveSurface( nSurface ) )
  {
    this->SendBroadcast( "LayerActorUpdated", this );
  }
}

int LayerSurface::GetActiveSurface()
{
  return ( m_surfaceSource ? m_surfaceSource->GetActiveSurface() : -1 );
}

void LayerSurface::SetActiveVector( int nVector )
{
  if ( m_surfaceSource && m_surfaceSource->SetActiveVector( nVector ) )
  {
    SetVisible( IsVisible() );  // refresh visibility
  }
}

int LayerSurface::GetActiveVector()
{
  return ( m_surfaceSource ? m_surfaceSource->GetActiveSurface() : -1 );
}

int LayerSurface::GetNumberOfVectorSets()
{
  return ( m_surfaceSource ? m_surfaceSource->GetNumberOfVectorSets() : 0 );
}

void LayerSurface::GetCurvatureRange( double* range )
{
  if ( m_surfaceSource && m_surfaceSource->GetMRIS() )
  {
    range[0] = m_surfaceSource->GetMRIS()->min_curv;
    range[1] = m_surfaceSource->GetMRIS()->max_curv;
  } 
}

bool LayerSurface::HasOverlay()
{
  return m_overlays.size() > 0;
}

void LayerSurface::SetActiveOverlay( int nOverlay )
{
  if ( nOverlay < (int)m_overlays.size() )
  {
    if ( m_nActiveOverlay < 0 && nOverlay >= 0 )
    {
      this->BlockListen( true );
      this->GetProperties()->SetCurvatureMap( LayerPropertiesSurface::CM_Binary );
      this->BlockListen( false );
    }
    m_nActiveOverlay = nOverlay;
    UpdateOverlay();
    this->SendBroadcast( "ActiveOverlayChanged", this );
  }
}

void LayerSurface::SetActiveOverlay( const char* name )
{
  for ( size_t i = 0; i < m_overlays.size(); i++ )
  {
    if ( strcmp( m_overlays[i]->GetName(), name ) == 0 )
    {
      SetActiveOverlay( i );
      return;
    }
  }
}

int LayerSurface::GetNumberOfOverlays()
{
  return m_overlays.size();
}
  
SurfaceOverlay* LayerSurface::GetOverlay( const char* name )
{
  for ( size_t i = 0; i < m_overlays.size(); i++ )
  {
    if ( strcmp( m_overlays[i]->GetName(), name ) == 0 )
    {
      return m_overlays[i];
    }
  }  
  
  return NULL;
}

int LayerSurface::GetActiveOverlayIndex()
{
  return m_nActiveOverlay;
}

SurfaceOverlay* LayerSurface::GetActiveOverlay()
{
  if ( m_nActiveOverlay >= 0 )
    return m_overlays[ m_nActiveOverlay ];
  else
    return NULL;
}

SurfaceOverlay* LayerSurface::GetOverlay( int n )
{
  if ( n >= 0 && n < (int)m_overlays.size() )
    return m_overlays[n];
  else
    return NULL;
}

void LayerSurface::UpdateOverlay( bool bAskRedraw )
{  
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mainActor->GetMapper() );
  vtkPolyData* polydata = mapper->GetInput();
  if ( m_nActiveOverlay >= 0 )
  {
    if ( mapper )
    {
      int nCount = polydata->GetPoints()->GetNumberOfPoints();
      this->BlockListen( true );
//      if ( mProperties->GetCurvatureMap() == LayerPropertiesSurface::CM_Threshold )
//        mProperties->SetCurvatureMap( LayerPropertiesSurface::CM_Binary );
      this->BlockListen( false );
      vtkSmartPointer<vtkUnsignedCharArray> array = vtkUnsignedCharArray::SafeDownCast( polydata->GetPointData()->GetArray( "Overlay" ) );
      if ( array.GetPointer() == NULL )
      { 
          array = vtkSmartPointer<vtkUnsignedCharArray>::New(); 
          array->SetNumberOfComponents( 4 );  
          array->SetNumberOfTuples( nCount );   
          array->SetName( "Overlay" );  
          polydata->GetPointData()->AddArray( array );
      }
      unsigned char* data = new unsigned char[ nCount*4 ];
      mProperties->GetCurvatureLUT()->MapScalarsThroughTable( polydata->GetPointData()->GetScalars("Curvature"), data, VTK_RGBA );
      GetActiveOverlay()->MapOverlay( data );
      for ( int i = 0; i < nCount; i++ )
      {
        array->SetTupleValue( i, data + i*4 );
      } 
      delete[] data;
      polydata->GetPointData()->SetActiveScalars( "Overlay" );
    }
  }
  else
  {
    polydata->GetPointData()->SetActiveScalars( "Curvature" );
  }
  if ( bAskRedraw )
    this->SendBroadcast( "LayerActorUpdated", this );
}
