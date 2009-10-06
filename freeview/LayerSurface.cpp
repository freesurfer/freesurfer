/**
 * @file  LayerSurface.cpp
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/10/06 21:46:47 $
 *    $Revision: 1.30 $
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
#include "vtkRGBAColorTransferFunction.h"
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
#include "SurfaceAnnotation.h"

LayerSurface::LayerSurface( LayerMRI* ref ) : Layer(),
    m_surfaceSource( NULL ),
    m_bResampleToRAS( true ),
    m_volumeRef( ref ),
    m_nActiveOverlay( -1 ),
    m_nActiveAnnotation( -1 )
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
  m_mainActor = vtkSmartPointer<vtkActor>::New();
  m_mainActor->GetProperty()->SetEdgeColor( 0.75, 0.75, 0.75 );
  mLowResFilter = vtkSmartPointer<vtkDecimatePro>::New();
  mLowResFilter->SetTargetReduction( 0.9 );
// mMediumResFilter = vtkSmartPointer<vtkDecimatePro>::New();
// mMediumResFilter->SetTargetReduction( 0.9 );

  m_vectorActor = vtkSmartPointer<vtkActor>::New();
  m_vectorActor->GetProperty()->SetColor( mProperties->GetVectorColor() );
  m_vectorActor->GetProperty()->SetPointSize( mProperties->GetVectorPointSize() );
  
  m_vertexActor = vtkSmartPointer<vtkActor>::New();
  m_vertexActor->GetProperty()->SetRepresentationToPoints();
  m_vertexActor->VisibilityOff();
  
  m_wireframeActor = vtkSmartPointer<vtkActor>::New();
  m_wireframeActor->VisibilityOff();
}

LayerSurface::~LayerSurface()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->Delete();
    m_sliceActor3D[i]->Delete();
  }

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


bool LayerSurface::LoadAnnotationFromFile( const char* filename )
{
  // create annotation
  SurfaceAnnotation* annot = new SurfaceAnnotation( this ); 
  bool ret = annot->LoadAnnotation( filename );
  if ( !ret )
  {
    delete annot;
    return false;
  }
  
  wxString fn = filename;
  fn = fn.substr( fn.find_last_of("/\\")+1 ).c_str();
  if ( fn.Right( 6 ) == ".annot" )
    fn = fn.Left( fn.Length()-6 );
  annot->SetName( fn.c_str() );
  
  m_annotations.push_back( annot );
  
  SetActiveAnnotation( m_annotations.size() - 1 );
  
  this->SendBroadcast( "LayerModified", this );
  this->SendBroadcast( "LayerAnnotationAdded", this );
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
  mapper->SetInput(  m_surfaceSource->GetVectorPolyData() );
  m_vectorActor->SetMapper( mapper );
//  mapper->Update();
  
  // vertex actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(  m_surfaceSource->GetVertexPolyData() );
  m_vertexActor->SetMapper( mapper );
  
  // wireframe actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(  m_surfaceSource->GetWireframePolyData() );
  m_wireframeActor->SetMapper( mapper );
  mapper->Update();

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
  m_wireframeActor->GetProperty()->SetColor( mProperties->GetBinaryColor() );
  m_vectorActor->GetProperty()->SetColor( mProperties->GetVectorColor() );
  if ( m_surfaceSource->IsCurvatureLoaded() )
  {
    if ( mProperties->GetCurvatureLUT() != m_mainActor->GetMapper()->GetLookupTable() )
    {
      m_mainActor->GetMapper()->SetLookupTable( mProperties->GetCurvatureLUT() );
      
    }
    
    if ( mProperties->GetMeshColorMap() == LayerPropertiesSurface::MC_Surface )
      m_wireframeActor->GetMapper()->SetLookupTable( mProperties->GetCurvatureLUT() );
    else if ( mProperties->GetMeshColorMap() == LayerPropertiesSurface::MC_Curvature )
      UpdateMeshRender();
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
  for ( int i = 0; i < 3; i++ )
    renderer->AddViewProp( m_sliceActor3D[i] );
  
  renderer->AddViewProp( m_mainActor );
  renderer->AddViewProp( m_vectorActor );
  renderer->AddViewProp( m_vertexActor );
  renderer->AddViewProp( m_wireframeActor );
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
  else if ( iMessage == "SurfaceRenderModeChanged" )
  {
    this->UpdateRenderMode();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "VertexRenderChanged" )
  {
    this->UpdateVertexRender();
    this->SendBroadcast( "LayerActorUpdated", this );
  }
  else if ( iMessage == "MeshRenderChanged" )
  {
    this->UpdateMeshRender();
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
  
  m_mainActor->SetVisibility( bVisible && mProperties->GetSurfaceRenderMode() != LayerPropertiesSurface::SM_Wireframe );
  m_wireframeActor->SetVisibility( bVisible && mProperties->GetSurfaceRenderMode() != LayerPropertiesSurface::SM_Surface );
  m_vertexActor->SetVisibility( bVisible && mProperties->GetShowVertices() );
  m_vectorActor->SetVisibility( ( bVisible && m_surfaceSource && m_surfaceSource->GetActiveVector() >= 0 )? 1 : 0 );

  this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerSurface::IsVisible()
{
  return m_sliceActor2D[0]->GetVisibility() > 0;
}

bool LayerSurface::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( m_sliceActor2D[i] == prop || m_sliceActor3D[i] == prop )
      return true;
  }
  return (m_mainActor.GetPointer() == prop);
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

double LayerSurface::GetCurvatureValue( int nVertex )
{
  return m_surfaceSource->GetCurvatureValue( nVertex );
}

bool LayerSurface::HasCurvature()
{
  return m_surfaceSource->IsCurvatureLoaded();
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
  vtkPolyData* polydataWireframe = vtkPolyDataMapper::SafeDownCast( m_wireframeActor->GetMapper() )->GetInput();
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
          polydataWireframe->GetPointData()->AddArray( array );
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
      if ( mProperties->GetMeshColorMap() == LayerPropertiesSurface::MC_Surface )
        polydataWireframe->GetPointData()->SetActiveScalars( "Overlay" );
    }
  }
  else if ( m_nActiveAnnotation >= 0 )
  {
    UpdateAnnotation( false );
  }
  else
  {
    polydata->GetPointData()->SetActiveScalars( "Curvature" );
    if ( mProperties->GetMeshColorMap() == LayerPropertiesSurface::MC_Surface )
      polydataWireframe->GetPointData()->SetActiveScalars( "Curvature" );
  }
  if ( bAskRedraw )
    this->SendBroadcast( "LayerActorUpdated", this );
}

void LayerSurface::UpdateRenderMode()
{
//  m_mainActor->GetProperty()->EdgeVisibilityOff();
//  m_mainActor->GetProperty()->BackfaceCullingOn();
  switch ( GetProperties()->GetSurfaceRenderMode() )
  {
    case LayerPropertiesSurface::SM_Surface:
      m_mainActor->VisibilityOn();
      m_wireframeActor->VisibilityOff();
      break;
    case LayerPropertiesSurface::SM_Wireframe:
      m_mainActor->VisibilityOff();
      m_wireframeActor->VisibilityOn();
      m_wireframeActor->GetProperty()->SetLineWidth( 1 );
      break;
    case LayerPropertiesSurface::SM_SurfaceAndWireframe:
      m_mainActor->VisibilityOn();
      m_wireframeActor->VisibilityOn();     
      m_wireframeActor->GetProperty()->SetLineWidth( 2 );
      break;
  }
}

void LayerSurface::SetActiveAnnotation( int n )
{
  if ( n < (int)m_annotations.size() )
  {
    m_nActiveAnnotation = n;
    UpdateAnnotation();
    this->SendBroadcast( "ActiveAnnotationChanged", this );
  }
}

void LayerSurface::SetActiveAnnotation( const char* name )
{
  for ( size_t i = 0; i < m_annotations.size(); i++ )
  {
    if ( strcmp( m_annotations[i]->GetName(), name ) == 0 )
    {
      SetActiveAnnotation( i );
      return;
    }
  }
}

int LayerSurface::GetNumberOfAnnotations()
{
  return m_annotations.size();
}
  
SurfaceAnnotation* LayerSurface::GetAnnotation( const char* name )
{
  for ( size_t i = 0; i < m_annotations.size(); i++ )
  {
    if ( strcmp( m_annotations[i]->GetName(), name ) == 0 )
    {
      return m_annotations[i];
    }
  }  
  
  return NULL;
}

int LayerSurface::GetActiveAnnotationIndex()
{
  return m_nActiveAnnotation;
}

SurfaceAnnotation* LayerSurface::GetActiveAnnotation()
{
  if ( m_nActiveAnnotation >= 0 )
    return m_annotations[ m_nActiveAnnotation ];
  else
    return NULL;
}

SurfaceAnnotation* LayerSurface::GetAnnotation( int n )
{
  if ( n >= 0 && n < (int)m_annotations.size() )
    return m_annotations[n];
  else
    return NULL;
}

void LayerSurface::UpdateAnnotation( bool bAskRedraw )
{  
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mainActor->GetMapper() );
  vtkPolyData* polydata = mapper->GetInput();
  vtkPolyDataMapper* mapperWireframe = vtkPolyDataMapper::SafeDownCast( m_wireframeActor->GetMapper() );
  vtkPolyData* polydataWireframe = mapperWireframe->GetInput();
  if ( m_nActiveAnnotation >= 0 )
  {
    if ( mapper )
    {
      int nCount = polydata->GetPoints()->GetNumberOfPoints();
      vtkSmartPointer<vtkIntArray> array = vtkIntArray::SafeDownCast( polydata->GetPointData()->GetArray( "Annotation" ) );
      if ( array.GetPointer() == NULL )
      { 
        array = vtkSmartPointer<vtkIntArray>::New();          
     //   array->SetNumberOfTuples( nCount ); 
        array->SetName( "Annotation" );  
        polydata->GetPointData()->AddArray( array );
        polydataWireframe->GetPointData()->AddArray( array );
      }

      array->SetArray( GetActiveAnnotation()->GetIndices(), nCount, 1 );
      polydata->GetPointData()->SetActiveScalars( "Annotation" );
      
      vtkSmartPointer<vtkFreesurferLookupTable> lut = vtkSmartPointer<vtkFreesurferLookupTable>::New();
      lut->BuildFromCTAB( GetActiveAnnotation()->GetColorTable(), false );  // do not clear zero
      mapper->SetLookupTable( lut );
      mapper->UseLookupTableScalarRangeOn();
      if ( mProperties->GetMeshColorMap() == LayerPropertiesSurface::MC_Surface )
      {
        polydataWireframe->GetPointData()->SetActiveScalars( "Annotation" );
        mapperWireframe->SetLookupTable( lut );
        mapperWireframe->UseLookupTableScalarRangeOn();
      }
    }
  }
  else
  {
    UpdateColorMap();
  }
  if ( bAskRedraw )
    this->SendBroadcast( "LayerActorUpdated", this );
}

void LayerSurface::UpdateVertexRender()
{
  m_vertexActor->SetVisibility( mProperties->GetShowVertices()? 1: 0 );
  m_vertexActor->GetProperty()->SetPointSize( mProperties->GetVertexPointSize() );
  m_vertexActor->GetProperty()->SetColor( mProperties->GetVertexColor() );
}

void LayerSurface::UpdateMeshRender()
{
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mainActor->GetMapper() );
  vtkPolyData* polydata = mapper->GetInput();
  vtkPolyDataMapper* mapperWireframe = vtkPolyDataMapper::SafeDownCast( m_wireframeActor->GetMapper() );
  vtkPolyData* polydataWireframe = mapperWireframe->GetInput();
  mapperWireframe->SetScalarVisibility( mProperties->GetMeshColorMap() != LayerPropertiesSurface::MC_Solid ? 1:0 );
  switch ( mProperties->GetMeshColorMap() )
  {
    case LayerPropertiesSurface::MC_Surface:
      polydataWireframe->GetPointData()->SetActiveScalars( polydata->GetPointData()->GetScalars()->GetName() );
      mapperWireframe->SetLookupTable( mapper->GetLookupTable() );
      break;
    case LayerPropertiesSurface::MC_Curvature:
      {
        // always display as threshold for curvature
        vtkSmartPointer<vtkRGBAColorTransferFunction> lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
        mProperties->BuildCurvatureLUT( lut, LayerPropertiesSurface::CM_Threshold );
        m_wireframeActor->GetMapper()->SetLookupTable( lut );
        polydataWireframe->GetPointData()->SetActiveScalars( "Curvature" );
      }
      break;
    case LayerPropertiesSurface::MC_Overlay:
      polydataWireframe->GetPointData()->SetActiveScalars( "Overlay" );
      break;
    case LayerPropertiesSurface::MC_Solid:
      m_wireframeActor->GetProperty()->SetColor( mProperties->GetMeshColor() );
      break;
    default:
      break;
  }
}
