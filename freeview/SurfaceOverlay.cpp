/**
 * @file  SurfaceOverlay.cxx
 * @brief Implementation for surface layer properties.
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
 *    $Date: 2011/03/12 00:28:53 $
 *    $Revision: 1.8 $
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


#include "SurfaceOverlay.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
#include "SurfaceOverlayProperty.h"
#include "FSSurface.h"
#include <QDebug>

SurfaceOverlay::SurfaceOverlay ( LayerSurface* surf ) :
    QObject(),
    m_fData( NULL ),
    m_dMaxValue(0),
    m_dMinValue(0),
    m_surface( surf ),
    m_bCorrelationData( false ),
    m_mriCorrelation(0),
    m_overlayPaired(0)
{
  InitializeData();  
  
  m_property =  new SurfaceOverlayProperty( this );
  connect( m_property, SIGNAL(ColorMapChanged()), surf, SLOT(UpdateOverlay()), Qt::UniqueConnection);
}

SurfaceOverlay::~SurfaceOverlay ()
{
  if ( m_fData )
    delete[] m_fData;

  if (m_overlayPaired)
  {
      m_overlayPaired->m_overlayPaired = 0;
  }
  else
  {
    delete m_property;
    if (m_mriCorrelation)
        MRIfree(&m_mriCorrelation);
  }
}

void SurfaceOverlay::InitializeData()
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    if ( m_fData )
      delete[] m_fData;
    m_nDataSize = mris->nvertices;
    m_fData = new float[ m_nDataSize ];
    if ( !m_fData )
      return;

    m_dMaxValue = m_dMinValue = mris->vertices[0].val;
    for ( int vno = 0; vno < m_nDataSize; vno++ )
    {
      m_fData[vno] = mris->vertices[vno].val;
      if ( m_dMaxValue < m_fData[vno] )
        m_dMaxValue = m_fData[vno];
      else if ( m_dMinValue > m_fData[vno] )
        m_dMinValue = m_fData[vno];
    }
  }
}

void SurfaceOverlay::CopyCorrelationData(SurfaceOverlay *overlay)
{
    if (!overlay->HasCorrelationData())
        return;
    delete m_property;
    m_property = overlay->m_property;
    connect( m_property, SIGNAL(ColorMapChanged()), m_surface, SLOT(UpdateOverlay()), Qt::UniqueConnection);
    m_mriCorrelation = overlay->m_mriCorrelation;
    m_overlayPaired = overlay;
    overlay->m_overlayPaired = this;
    m_bCorrelationData = true;
}

bool SurfaceOverlay::LoadCorrelationData( const QString& filename )
{
  MRI* mri = ::MRIreadHeader( filename.toAscii().data(), -1 );
  if ( mri == NULL )
  {
    cerr << "MRIread failed: unable to read from " << qPrintable(filename) << "\n";
    return false;
  }
  if ( mri->width != m_nDataSize*2 || (mri->height != 1 && mri->height != m_nDataSize*2) ||
       (mri->nframes != 1 && mri->nframes != m_nDataSize*2))
  {
    cerr << "Correlation data does not match with surface\n";
    MRIfree( &mri );
    return false;
  }
  MRIfree( &mri );
  mri = ::MRIread( filename.toAscii().data() );      // long process
  if ( mri == NULL )
  {
    cerr << "MRIread failed: Unable to read from " << qPrintable(filename) << "\n";
    return false;
  }
  m_mriCorrelation = mri;
  m_bCorrelationData = true;
  m_bCorrelationDataReady = false;
  
  return true;
}

void SurfaceOverlay::UpdateCorrelationAtVertex( int nVertex, int nHemisphere )
{
  if ( nHemisphere == -1)
      nHemisphere = m_surface->GetHemisphere();
  int nVertexOffset = nHemisphere * m_nDataSize;
  int nDataOffset = m_surface->GetHemisphere() * m_nDataSize;
  double old_range = m_dMaxValue - m_dMinValue;
  if (m_mriCorrelation->height > 1)
    m_dMaxValue = m_dMinValue =
                  MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, nDataOffset, 0, 0 );
  else
    m_dMaxValue = m_dMinValue =
                    MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, 0, 0, nDataOffset );
  for ( int i = 0; i < m_nDataSize; i++ )
  {
    if (m_mriCorrelation->height > 1)
        m_fData[i] = MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, i + nDataOffset, 0, 0 );
    else
        m_fData[i] = MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, 0, 0, i + nDataOffset );
    if ( m_dMaxValue < m_fData[i] )
      m_dMaxValue = m_fData[i];
    else if ( m_dMinValue > m_fData[i] )
      m_dMinValue = m_fData[i];
  }
  m_bCorrelationDataReady = true;
  if (old_range <= 0)
      m_property->Reset();

  if (m_overlayPaired && nHemisphere == m_surface->GetHemisphere())
  {
      m_overlayPaired->blockSignals(true);
      m_overlayPaired->UpdateCorrelationAtVertex(nVertex, nHemisphere);
      m_overlayPaired->blockSignals(false);
  }

  m_surface->UpdateOverlay(true);
  emit DataUpdated();
}

QString SurfaceOverlay::GetName()
{
  return m_strName;
}

void SurfaceOverlay::SetName( const QString& name )
{
  m_strName = name;
}

void SurfaceOverlay::MapOverlay( unsigned char* colordata )
{
  if ( !m_bCorrelationData || m_bCorrelationDataReady )
    m_property->MapOverlayColor( m_fData, colordata, m_nDataSize );
}

double SurfaceOverlay::GetDataAtVertex( int nVertex )
{
  return m_fData[nVertex];
}


