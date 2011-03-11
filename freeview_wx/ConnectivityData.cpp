/**
 * @file  ConnectivityData.cpp
 * @brief Holder for connectivity data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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

#include "ConnectivityData.h"
#include "LayerSurface.h"
#include "LayerPropertiesSurface.h"
#include "SurfaceAnnotation.h"
#include <wx/wx.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkSplineFilter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkAppendPolyData.h>
#include <vtkRGBAColorTransferFunction.h>
#include <vtkRenderer.h>
#include <wx/file.h>
#include <wx/filename.h>

using namespace std;

ConnectivityData::ConnectivityData() : 
    Listener( "ConnectivityData" ), 
    Broadcaster( "ConnectivityData" ),
    m_bVaryRadius( true ),
    m_dRadius( 0.5 ),
    m_dRadiusMax( 6 ),
    m_dThresholdMin( 0 ),
    m_dThresholdMax( 1 ),
    m_nDisplayMode( 0 ),
    m_bIncrementalDisplay( false ),
    m_MRI( NULL ),
    m_lut( NULL ),
    m_surfLeft( NULL ),
    m_surfRight( NULL ),
    m_dMatConnectivity( NULL ),
    m_dMinConnectivity( 0 ),
    m_dMaxConnectivity( 1 ),
    m_nSeedLabels( NULL )
{
  m_actorStartPoint = vtkSmartPointer<vtkActor>::New();
  m_actorTube = vtkSmartPointer<vtkActor>::New();
}

ConnectivityData::~ConnectivityData()
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );
  if ( m_lut )
    ::CTABfree( &m_lut );
  
  if ( m_dMatConnectivity )
  {
    for ( int i = 0; i < m_nMatSize; i++ )
      delete[] m_dMatConnectivity[i];
    delete[] m_dMatConnectivity;
  }
  if ( m_nSeedLabels )
    delete[] m_nSeedLabels;
}
    
bool ConnectivityData::IsValid()
{
  return m_MRI && m_lut && m_surfLeft && m_surfRight;
}

void ConnectivityData::AppendProps( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorStartPoint );
  renderer->AddViewProp( m_actorTube );
}

bool ConnectivityData::Load( const char* fn, 
                             const char* lut_fn, 
                             LayerSurface* surf1, 
                             LayerSurface* surf2 )
{
  if ( m_MRI )
    ::MRIfree( &m_MRI );
  
  m_MRI = ::MRIread( fn );
  if ( !m_MRI )
  {
    cerr << "Can not load connectivity file " << fn << endl;
    return false;
  }
  
  m_name = wxFileName( fn ).GetFullName();
  
  if ( surf1->GetActiveAnnotation() < 0 || surf2->GetActiveAnnotation() < 0 )
  {
    cerr << "No surface annotation was found. Can not display connectivity." << endl;
    return false;
  }  
  
  if ( ( surf1->GetHemisphere() == 0 && surf2->GetHemisphere() == 1 ) ||
       ( surf1->GetHemisphere() == 1 && surf2->GetHemisphere() == 0 ) )
  { 
    m_surfLeft = surf1->GetHemisphere() ? surf2 : surf1;
    m_surfRight = surf1->GetHemisphere() ? surf1 : surf2;
  }
  else
  {
    cerr << "One of the hemispheres is missing!" << endl;
    return false;
  }

  m_surfLeft->AddListener( this );
  m_surfRight->AddListener( this );
  
  if ( m_lut )
    ::CTABfree( &m_lut );
  
  m_lut = CTABreadASCII( lut_fn );
  if ( !m_lut )
  {
    m_lut = CTABdeepCopy( surf1->GetAnnotation(0)->GetColorTable() );
    if ( !m_lut )
    { 
      cerr << "Can not load lookup table file " << lut_fn << endl;
      return false;
    }
  }
    
  m_dMinConnectivity = m_dMaxConnectivity = MRIgetVoxVal( m_MRI, 0, 0, 0, 0 );
  
  for ( int i = 0; i < m_MRI->width; i++ )
  {
    for ( int m = 0; m < m_MRI->height; m++ )
    for ( int k = 0; k < m_MRI->depth; k++ )
    for ( int j = 0; j < m_MRI->nframes; j++ )
    {
      double val = MRIgetVoxVal( m_MRI, i, m, k, j );
      if ( val > m_dMaxConnectivity )
        m_dMaxConnectivity = val;
      else if ( val < m_dMinConnectivity )
        m_dMinConnectivity = val;
    }
  }
  m_dThresholdMax = m_dMaxConnectivity;
  
  cout << "Connectivity value range: " << m_dMinConnectivity << ", " << m_dMaxConnectivity << endl;
  
  // build connectivity data between annotations
  if ( m_dMatConnectivity )
  {
    for ( int i = 0; i < m_nMatSize; i++ )
      delete[] m_dMatConnectivity[i];
    delete[] m_dMatConnectivity;
  }
  m_nMatSize = m_surfLeft->GetAnnotation( 0 )->GetNumberOfAnnotations() + m_surfRight->GetAnnotation( 0 )->GetNumberOfAnnotations();  

  m_dMatConnectivity = new double*[m_nMatSize];
  for ( int i = 0; i < m_nMatSize; i++ )
    m_dMatConnectivity[i] = new double[m_nMatSize];
  
  int nLeftSize = m_surfLeft->GetAnnotation( 0 )->GetNumberOfAnnotations();  
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    for ( int j = i+1; j < m_nMatSize; j++ )
    {
      int ni, nj;
      if ( i < nLeftSize )
        ni = GetColorTableIndex( m_surfLeft, i );
      else
        ni = GetColorTableIndex( m_surfRight, i-nLeftSize );
      if ( j < nLeftSize )
        nj = GetColorTableIndex( m_surfLeft, j );
      else
        nj = GetColorTableIndex( m_surfRight, j-nLeftSize );
          
      if ( ni >= 0 && nj >= 0 )
      {
        if ( m_MRI->depth <= 1 && m_MRI->nframes <= 1 )
          m_dMatConnectivity[i][j] = m_dMatConnectivity[j][i] = MRIgetVoxVal( m_MRI, ni, nj, 0, 0 );
        else
          m_dMatConnectivity[i][j] = m_dMatConnectivity[j][i] = MRIgetVoxVal( m_MRI, ni, 0, 0, nj );
      }
      else
      {
        m_dMatConnectivity[i][j] = m_dMatConnectivity[j][i] = m_dMinConnectivity - 0.01;  // make sure it's smaller than the minimum connectivity
      }
    }
  }
  
  if ( m_nSeedLabels )
    delete[] m_nSeedLabels;
  m_nSeedLabels = new int[m_nMatSize];
  memset( m_nSeedLabels, 0, sizeof(int)*m_nMatSize );
  
  // build lut table
  m_vtkLut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    int rgb[3];
    if ( i < nLeftSize )
      m_surfLeft->GetAnnotation(0)->GetAnnotationColorAtIndex( i, rgb );
    else
      m_surfRight->GetAnnotation(0)->GetAnnotationColorAtIndex( i-nLeftSize, rgb );
    m_vtkLut->AddRGBAPoint( i, rgb[0]/255.0, rgb[1]/255.0, rgb[2]/255.0, 1 );
  }
  m_vtkLut->Build();
  
  return true;
}

void ConnectivityData::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerRemoved" )
  {
    if ( iData == m_surfLeft || iData == m_surfRight )
    {
      m_surfLeft = NULL;
      m_surfRight = NULL;
    }
  }
}

void ConnectivityData::BuildConnectivityActors()
{
  if ( !IsValid() )
    return;
  
  int nLeftSize = m_surfLeft->GetAnnotation( 0 )->GetNumberOfAnnotations();
  if ( m_nDisplayMode == DM_All )
  {
    for ( int i = 0; i < m_nMatSize; i++ )
      m_nSeedLabels[i] = 1;
  }
  else 
  {
    if ( !m_bIncrementalDisplay )
      memset( m_nSeedLabels, 0, sizeof(int)*m_nMatSize );
    
    double* RASPos = m_surfLeft->GetSlicePosition();
    int index_annot = -1;
    int nVertex = m_surfLeft->GetVertexIndexAtTarget( RASPos, NULL );
    if ( nVertex >= 0 )
    {
      index_annot = m_surfLeft->GetAnnotation( 0 )->GetIndexAtVertex( nVertex );
    }
    else if ( ( nVertex = m_surfRight->GetVertexIndexAtTarget( RASPos, NULL ) ) >= 0 )
    {
      index_annot = m_surfRight->GetAnnotation( 0 )->GetIndexAtVertex( nVertex ) + nLeftSize;
    }
    
    if ( index_annot >= 0 )
    {
      m_nSeedLabels[index_annot] = 1;
      if ( m_nDisplayMode == DM_Both )
        m_nSeedLabels[(index_annot<nLeftSize?index_annot+nLeftSize:index_annot-nLeftSize)] = 1;
    }
    else 
    {
      if ( m_nLastDisplayMode == DM_All )
      {
        // clear actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
        m_actorTube->SetMapper( mapper );
        this->SendBroadcast( "ConnectivityActorUpdated", this );
      }
      return;
    }
  }
  
  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  int nAppendCount = 0;  
  double dPos1[3], dPos2[3];
  short* nMask = new short[m_nMatSize*m_nMatSize];
  memset( nMask, 0, sizeof(short)*m_nMatSize*m_nMatSize );
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    if ( m_nSeedLabels[i] )
    {
      for ( int j = 0; j < m_nMatSize; j++ )
      {
        if ( i != j && nMask[i*m_nMatSize+j] == 0 &&
             m_dMatConnectivity[i][j] >= m_dThresholdMin &&
             m_dMatConnectivity[i][j] <= m_dThresholdMax )
        {
          if ( i < nLeftSize )
            m_surfLeft->GetAnnotation(0)->GetAnnotationPoint( i, dPos1 );
          else
            m_surfRight->GetAnnotation(0)->GetAnnotationPoint( i-nLeftSize, dPos1 );
          if ( j < nLeftSize )
            m_surfLeft->GetAnnotation(0)->GetAnnotationPoint( j, dPos2 );
          else
            m_surfRight->GetAnnotation(0)->GetAnnotationPoint( j-nLeftSize, dPos2 );
          vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
          BuildWireActor( i, dPos1, dPos2, m_dMatConnectivity[i][j], polydata );
          append->AddInput( polydata );
          nAppendCount++;
          nMask[i*m_nMatSize+j] = nMask[j*m_nMatSize+i] = 1;
        }
      }
    }    
  }
  delete[] nMask;

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
  if ( nAppendCount > 0 )
  {
    vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
    spline->SetInputConnection( append->GetOutputPort() );
    spline->SetSubdivideToLength();
    spline->SetLength( 3 );
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInputConnection( spline->GetOutputPort() );
    tube->SetRadius( m_dRadius );
    tube->SetRadiusFactor( m_dRadiusMax / m_dRadius );
    tube->SetNumberOfSides( 6 );
    if ( m_bVaryRadius )
      tube->SetVaryRadiusToVaryRadiusByScalar();
    
    spline->GetOutput()->GetPointData()->SetActiveScalars( "connectivity" );
    tube->Update();
    mapper->SetInputConnection( tube->GetOutputPort() );
    mapper->SetLookupTable( m_vtkLut );
    mapper->ScalarVisibilityOn();
    mapper->GetInput()->GetPointData()->SetActiveScalars( "annot_index" );
  }
  m_actorTube->SetMapper( mapper );
  m_nLastDisplayMode = m_nDisplayMode;
  
  this->SendBroadcast( "ConnectivityActorUpdated", this );
}

void ConnectivityData::BuildWireActor( int seed_index, double* dStartPos, double* dPos, double val, vtkPolyData* polydata )
{
  double mid_pt[3];
  for ( int i = 0; i < 3; i++ )
  {
    mid_pt[i] = ( dStartPos[i] + dPos[i] ) / 2;
  }
  mid_pt[2] += sqrt( vtkMath::Distance2BetweenPoints( dStartPos, dPos ) ) / 2;
  
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  pts->InsertNextPoint( dStartPos );
  pts->InsertNextPoint( mid_pt );
  pts->InsertNextPoint( dPos );
  lines->InsertNextCell( 3 );
  lines->InsertCellPoint( 0 );
  lines->InsertCellPoint( 1 );
  lines->InsertCellPoint( 2 );
  scalars->SetNumberOfValues( 3 );
  scalars->SetValue( 0, val );
  scalars->SetValue( 1, val );
  scalars->SetValue( 2, val );
  scalars->SetName( "connectivity" );
  polydata->SetPoints( pts );
  polydata->SetLines( lines );
  polydata->GetPointData()->SetScalars( scalars );
  vtkSmartPointer<vtkFloatArray> scalars2 = vtkSmartPointer<vtkFloatArray>::New();
  scalars2->SetNumberOfValues( 3 );
  scalars2->SetValue( 0, seed_index );
  scalars2->SetValue( 1, seed_index );
  scalars2->SetValue( 2, seed_index );
  scalars2->SetName( "annot_index" );
  polydata->GetPointData()->AddArray( scalars2 );
}


bool ConnectivityData::HasAnySeeds()
{
  if ( !IsValid() )
    return false;
  
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    if ( m_nSeedLabels[i] )
      return true;
  }
  return false;
}

bool ConnectivityData::Export( const char* filename )
{
  wxFile file;
  if ( !file.Open( filename, wxFile::write ) )
  {
    cerr << "Failed to open file for write: " << filename << endl;
    return false;
  }
  
  file.Write( "graph connectivity {\n" );
  
  // write connection info to a string buffer first
  short* nMask = new short[m_nMatSize*m_nMatSize];
  memset( nMask, 0, sizeof(short)*m_nMatSize*m_nMatSize );
  short* nConns = new short[m_nMatSize];
  memset( nConns, 0, sizeof(short)*m_nMatSize );
  wxString strg_buffer = "";
  
  //  m_nMatSize /= 2;  /***** hack */
  
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    if ( m_nSeedLabels[i] )
    {
      for ( int j = 0; j < m_nMatSize; j++ )
      {
        if ( i != j && nMask[i*m_nMatSize+j] == 0 &&
             m_dMatConnectivity[i][j] >= m_dThresholdMin &&
             m_dMatConnectivity[i][j] <= m_dThresholdMax )
        {
          strg_buffer += wxString::Format( "%d -- %d [label=\"%.2f\", len=1.0];\n", i, j, m_dMatConnectivity[i][j] );
//          strg_buffer += wxString::Format( "%d -- %d [len=1.0];\n", i, j );
          nMask[i*m_nMatSize+j] = nMask[j*m_nMatSize+i] = 1;
          nConns[i] = nConns[j] = 1;
        }
      }
    }
  }
  
  // write node properties (skip non-connection nodes)
  int nLeftSize = m_surfLeft->GetAnnotation( 0 )->GetNumberOfAnnotations();
  for ( int i = 0; i < m_nMatSize; i++ )
  {
    if ( nConns[i] )
    {
      wxString name;
      int rgb[3] = { 255, 255, 255 };
      if ( i < nLeftSize )
      {
        name = wxString("lh-") + m_surfLeft->GetAnnotation(0)->GetAnnotationNameAtIndex( i ).c_str();
        m_surfLeft->GetAnnotation(0)->GetAnnotationColorAtIndex( i, rgb );
      }
      else
      {
        name = wxString("rh-") + m_surfRight->GetAnnotation(0)->GetAnnotationNameAtIndex( i-nLeftSize ).c_str();
        m_surfRight->GetAnnotation(0)->GetAnnotationColorAtIndex( i-nLeftSize, rgb );
      }
//      file.Write( wxString::Format( "node [label=\"%s\", style=\"filled\", color=\"#%02x%02x%02x\"]; %d;\n", name.c_str(), rgb[0], rgb[1], rgb[2], i ) );
      file.Write( wxString::Format( "node [style=\"filled\", color=\"#%02x%02x%02x\"]; %d;\n", rgb[0], rgb[1], rgb[2], i ) );
    }
  }
  
  // write buffered connection data
  file.Write( strg_buffer );
  file.Write( wxString::Format( "label=\"%s (%f, %f)\"\n;", m_name.c_str(), m_dThresholdMin, m_dThresholdMax ) );
  file.Write( "}\n" );
  file.Flush();
  file.Close();
  delete[] nMask;
  delete[] nConns;
  
  //  m_nMatSize *= 2;  /***** hack back*/
  
  return true;
}

int ConnectivityData::GetColorTableIndex( LayerSurface* surf, int annotation_index )
{
  SurfaceAnnotation* annot = surf->GetAnnotation( 0 );
  wxString name = _("lh");
  if ( surf->GetHemisphere() > 0 )
    name = _("rh");
  
  name = name + "-" + annot->GetAnnotationNameAtIndex( annotation_index ).c_str();
  
  int nEntry;
  if ( CTABfindEntryByName( m_lut, name.c_str(), &nEntry ) == 0 && nEntry >= 0 )
    return nEntry;
  
  name = wxString(_("ctx-") ) + name;
  if ( CTABfindEntryByName( m_lut, name.c_str(), &nEntry ) == 0 && nEntry >= 0 )
    return nEntry;
  
  name = annot->GetAnnotationNameAtIndex( annotation_index ).c_str();
  if ( CTABfindEntryByName( m_lut, name.c_str(), &nEntry ) == 0 && nEntry >= 0 )
    return nEntry;
  
  return annotation_index;
}

void ConnectivityData::SetVaryRadius( bool bVary )
{
  if ( m_bVaryRadius != bVary )
  {
    m_bVaryRadius = bVary;
    BuildConnectivityActors();
  }
}

void ConnectivityData::SetRadius( double dRadius )
{
  if ( dRadius != m_dRadius )
  {
    if ( !m_bVaryRadius && dRadius > m_dRadiusMax )
      m_dRadiusMax = dRadius;
    
    SetRadiusRange( dRadius, m_dRadiusMax );
  }
}
    
void ConnectivityData::SetRadiusMax( double dRadius )
{
  if ( dRadius != m_dRadiusMax )
  {
    SetRadiusRange( m_dRadius, dRadius );
  }
}

void ConnectivityData::SetRadiusRange( double dMin, double dMax )
{
  m_dRadius = dMin;
  m_dRadiusMax = dMax;
  
  BuildConnectivityActors();
}
    
void ConnectivityData::SetThresholdMin  ( double dMin )
{
  SetThresholdRange( dMin, m_dThresholdMax );
}
    
void ConnectivityData::SetThresholdMax  ( double dMax )
{
  SetThresholdRange( m_dThresholdMin, dMax );
}

void ConnectivityData::SetThresholdRange( double dMin, double dMax )
{
  if ( m_dThresholdMin != dMin || m_dThresholdMax != dMax )
  {
    m_dThresholdMin = dMin;
    m_dThresholdMax = dMax;
    
    BuildConnectivityActors();
  }
}

void ConnectivityData::SetDisplayMode ( int nMode )
{
  if ( m_nDisplayMode != nMode )
  {
    m_nDisplayMode = nMode;
    BuildConnectivityActors();
  }
}

void ConnectivityData::SetIncrementalDisplay( bool bIncre )
{
  if ( m_bIncrementalDisplay != bIncre )
  {
    m_bIncrementalDisplay = bIncre;
    BuildConnectivityActors();
  }
}

