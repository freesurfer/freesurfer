/**
 * @file  SurfaceAnnotation.cxx
 * @brief Implementation for surface annotation.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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


#include <assert.h>
#include "SurfaceAnnotation.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
//#include "SurfaceAnnotationProperties.h"
#include "FSSurface.h"
#include <wx/filename.h>

SurfaceAnnotation::SurfaceAnnotation ( LayerSurface* surf ) :
    Broadcaster( "SurfaceAnnotation" ),
    Listener( "SurfaceAnnotation" ),
    m_nIndices( NULL ),
    m_nCenterVertices( NULL ),
    m_lut( NULL ),
    m_surface( surf )
{
}

SurfaceAnnotation::~SurfaceAnnotation ()
{
  Reset();
}

void SurfaceAnnotation::Reset()
{
  if ( m_nIndices )
    delete[] m_nIndices;
  
  if ( m_nCenterVertices )
    delete[] m_nCenterVertices;
   
  m_nIndices = NULL;
  m_lut = NULL;
  m_nCenterVertices = NULL;
}

void SurfaceAnnotation::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->SendBroadcast( "AnnotationChanged", this );
  }
}

bool SurfaceAnnotation::LoadAnnotation( const char* fn )
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    m_nIndices = NULL;

    if ( MRISreadAnnotation( mris, fn ) != 0 )
    {
      cerr << "Could not load annotation from file " << fn << endl;
      return false;
    }
    else
    {
      Reset();
      m_lut = mris->ct;
      m_nIndexSize = mris->nvertices;
      m_nIndices = new int[m_nIndexSize]; 
      CTABgetNumberOfValidEntries( m_lut, &m_nAnnotations ); 
      m_nCenterVertices = new int[m_nAnnotations]; 
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        m_nCenterVertices[i] = -1;
      }      
      
      // convert annotations to lookup table indices
      // also find the "center vertex"
      double** pts = new double*[m_nAnnotations];
      int* vcount = new int[m_nAnnotations];
      memset( vcount, 0, sizeof( int )*m_nAnnotations );
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        pts[i] = new double[3];
        memset( pts[i], 0, sizeof( double ) * 3 );
      }
      
      int n;
      for ( int i = 0; i < m_nIndexSize; i++ )
      {
        if ( CTABfindAnnotation( m_lut, mris->vertices[i].annotation, &n ) == 0 )
        {
          m_nIndices[i] = n;
          
          if ( mris->vertices[i].ripflag == 0 &&
               m_nIndices[i] != -1 )
          {
            vcount[n]++;
            pts[n][0] += mris->vertices[i].x;
            pts[n][1] += mris->vertices[i].y;
            pts[n][2] += mris->vertices[i].z;
            m_nCenterVertices[n] = i;
          }
        }
        else
          m_nIndices[i] = -1;
      }  
      
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        if ( vcount[i] > 0 )
        {
          pts[i][0] /= vcount[i];
          pts[i][1] /= vcount[i];
          pts[i][2] /= vcount[i];
          
          int nVertex = m_surface->GetSourceSurface()->FindVertexAtSurfaceRAS( pts[i], NULL );
          if ( nVertex >=0 && m_nIndices[nVertex] == i )
            m_nCenterVertices[i] = nVertex;
        }
      }
      
      delete[] vcount;
      for ( int i = 0; i < m_nAnnotations; i++ )
        delete[] pts[i];
      delete[] pts;
      
      return true;
    }
  }
  return false;
}

void SurfaceAnnotation::GetAnnotationPoint( int nIndex, double* pt_out )
{
  m_surface->GetTargetAtVertex( m_nCenterVertices[nIndex], pt_out );
}

const char* SurfaceAnnotation::GetName()
{
  return m_strName.c_str();
}

void SurfaceAnnotation::SetName( const char* name )
{
  m_strName = name;
}

int SurfaceAnnotation::GetIndexAtVertex( int nVertex )
{
  return m_nIndices[nVertex];
}

std::string SurfaceAnnotation::GetAnnotationNameAtVertex( int nVertex )
{
  return GetAnnotationNameAtIndex( m_nIndices[nVertex] );
}

std::string SurfaceAnnotation::GetAnnotationNameAtIndex( int nIndex )
{
  char name[128];
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
    CTABisEntryValid( m_lut, nIndex, &nValid );
  if ( nValid && CTABcopyName( m_lut, nIndex, name, 128 ) == 0 )
  {
    return name;
  }
  
  return ""; 
}

void SurfaceAnnotation::GetAnnotationColorAtIndex( int nIndex, int* rgb )
{
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
    CTABisEntryValid( m_lut, nIndex, &nValid );
  if ( nValid )
    CTABrgbAtIndexi( m_lut, nIndex, rgb, rgb+1, rgb+2 );
}
