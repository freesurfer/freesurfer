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
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/08/21 19:57:52 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2007-2009,
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


#include <assert.h>
#include "SurfaceAnnotation.h"
#include "vtkLookupTable.h"
#include "vtkRGBATransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
//#include "SurfaceAnnotationProperties.h"
#include "FSSurface.h"
#include <wx/filename.h>

using namespace std;

SurfaceAnnotation::SurfaceAnnotation ( LayerSurface* surf ) :
    Broadcaster( "SurfaceAnnotation" ),
    Listener( "SurfaceAnnotation" ),
    m_nIndices( NULL ),
    m_surface( surf ),
    m_lut( NULL )
{
}

SurfaceAnnotation::~SurfaceAnnotation ()
{
  if ( m_nIndices )
    delete[] m_nIndices;
  if ( m_lut )
    CTABfree( &m_lut );
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
      if ( m_lut )
        CTABfree( &m_lut );
      m_lut = mris->ct;
      
      if ( m_nIndices )
        delete( m_nIndices );
      m_nIndexSize = mris->nvertices;
      m_nIndices = new int[m_nIndexSize];  
      int n;
      for ( int i = 0; i < m_nIndexSize; i++ )
      {
        if ( CTABfindAnnotation( m_lut, mris->vertices[i].annotation, &n ) == 0 )
          m_nIndices[i] = n;
      }      
      
      return true;
    }
  }
  return false;
}

const char* SurfaceAnnotation::GetName()
{
  return m_strName.c_str();
}

void SurfaceAnnotation::SetName( const char* name )
{
  m_strName = name;
}

void SurfaceAnnotation::MapOverlay( unsigned char* colordata )
{
//  m_properties->MapOverlayColor( colordata, m_nDataSize );
}

int SurfaceAnnotation::GetIndexAtVertex( int nVertex )
{
  return m_nIndices[nVertex];
}

std::string SurfaceAnnotation::GetAnnotationNameAtVertex( int nVertex )
{
  char name[128];
  int nValid = 0;
  int nTotalCount = 0;
  int nIndex = m_nIndices[nVertex];
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
    CTABisEntryValid( m_lut, nIndex, &nValid );
  if ( nValid && CTABcopyName( m_lut, nIndex, name, 128 ) == 0 )
  {
    return name;
  }
  
  return ""; 
}


