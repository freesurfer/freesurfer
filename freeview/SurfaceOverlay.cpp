/**
 * @file  SurfaceOverlay.cpp
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
 *    $Author: rpwang $
 *    $Date: 2012/10/31 20:10:12 $
 *    $Revision: 1.15 $
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
#include "ProgressCallback.h"

extern "C"
{
#include "utils.h"
}

SurfaceOverlay::SurfaceOverlay ( LayerSurface* surf ) :
  QObject(),
  m_fData( NULL ),
  m_fDataRaw( NULL ),
  m_fDataOriginal( NULL ),
  m_dMaxValue(0),
  m_dMinValue(0),
  m_surface( surf ),
  m_bCorrelationData( false ),
  m_mriCorrelation(NULL),
  m_overlayPaired(NULL),
  m_nActiveFrame(0),
  m_nNumOfFrames(1)
{
  InitializeData();

  m_property =  new SurfaceOverlayProperty( this );
  connect( m_property, SIGNAL(ColorMapChanged()), surf, SLOT(UpdateOverlay()), Qt::UniqueConnection);
  connect( m_property, SIGNAL(SmoothChanged()), this, SLOT(UpdateSmooth()), Qt::UniqueConnection);
}

SurfaceOverlay::~SurfaceOverlay ()
{
  if ( m_fDataRaw )
    delete[] m_fDataRaw;

  if (m_fDataOriginal)
    delete[] m_fDataOriginal;

  if (m_overlayPaired)
  {
    m_overlayPaired->m_overlayPaired = 0;
  }
  else
  {
    delete m_property;
    if (m_mriCorrelation)
    {
      MRIfree(&m_mriCorrelation);
    }
  }
}

void SurfaceOverlay::InitializeData()
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    if ( m_fDataRaw )
      delete[] m_fDataRaw;

    m_nDataSize = mris->nvertices;
    m_fData = m_fDataRaw = new float[ m_nDataSize ];
    if ( !m_fData )
    {
      return;
    }

    if ( m_fDataOriginal )
      delete[] m_fDataOriginal;

    m_fDataOriginal = new float[ m_nDataSize ];
    if ( !m_fDataOriginal )
    {
      return;
    }

    m_dMaxValue = m_dMinValue = mris->vertices[0].val;
    for ( int vno = 0; vno < m_nDataSize; vno++ )
    {
      m_fData[vno] = mris->vertices[vno].val;
      if ( m_dMaxValue < m_fData[vno] )
      {
        m_dMaxValue = m_fData[vno];
      }
      else if ( m_dMinValue > m_fData[vno] )
      {
        m_dMinValue = m_fData[vno];
      }
    }
    memcpy(m_fDataOriginal, m_fData, sizeof(float)*m_nDataSize);
  }
}

void SurfaceOverlay::InitializeData(float *data_buffer_in, int nvertices, int nframes)
{
  if ( m_surface )
  {
    if ( m_fDataRaw )
      delete[] m_fDataRaw;

    m_nDataSize = nvertices;
    m_nNumOfFrames = nframes;
    m_fData = m_fDataRaw = data_buffer_in;
    if ( !m_fData )
    {
      return;
    }

    if ( m_fDataOriginal )
      delete[] m_fDataOriginal;

    m_fDataOriginal = new float[ m_nDataSize ];
    if ( !m_fDataOriginal )
    {
      return;
    }

    SetActiveFrame(0);
    GetProperty()->Reset();
  }
}

void SurfaceOverlay::CopyCorrelationData(SurfaceOverlay *overlay)
{
  if (!overlay->HasCorrelationData())
  {
    return;
  }
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
  ::SetProgressCallback(ProgressCallback, 0, 100);
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
  {
    nHemisphere = m_surface->GetHemisphere();
  }
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
    {
      m_fData[i] = MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, i + nDataOffset, 0, 0 );
    }
    else
    {
      m_fData[i] = MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, 0, 0, i + nDataOffset );
    }
    if ( m_dMaxValue < m_fData[i] )
    {
      m_dMaxValue = m_fData[i];
    }
    else if ( m_dMinValue > m_fData[i] )
    {
      m_dMinValue = m_fData[i];
    }
  }
  memcpy(m_fDataOriginal, m_fData, sizeof(float)*m_nDataSize);
  if (GetProperty()->GetSmooth())
  {
    SmoothData();
  }

  m_bCorrelationDataReady = true;
  if (old_range <= 0)
  {
    m_property->Reset();
  }

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
  {
    m_property->MapOverlayColor( m_fData, colordata, m_nDataSize );
  }
}

double SurfaceOverlay::GetDataAtVertex( int nVertex )
{
  return m_fData[nVertex];
}

void SurfaceOverlay::UpdateSmooth(bool trigger_paired)
{
  if (GetProperty()->GetSmooth())
  {
    SmoothData();
  }
  else
  {
    memcpy(m_fData, m_fDataOriginal, sizeof(float)*m_nDataSize);
  }
  m_surface->UpdateOverlay(true);
  emit DataUpdated();

  if (trigger_paired && m_overlayPaired)
    m_overlayPaired->UpdateSmooth(false);
}

void SurfaceOverlay::SmoothData(int nSteps_in, float *data_out)
{
  MRI* mri = MRIallocSequence( m_nDataSize,
                               1,
                               1,
                               MRI_FLOAT, 1);
  if (!mri)
  {
    cerr << "Can not allocate mri\n";
    return;
  }
  memcpy(&MRIFseq_vox( mri, 0, 0, 0, 0 ), m_fDataOriginal, sizeof(float)*m_nDataSize);
  int nSteps = nSteps_in;
  if (nSteps < 1)
    nSteps = GetProperty()->GetSmoothSteps();
  MRI* mri_smoothed = ::MRISsmoothMRI(m_surface->GetSourceSurface()->GetMRIS(), mri,
                                      nSteps, NULL, NULL);
  if (mri_smoothed)
  {
    if (data_out)
      memcpy(data_out, &MRIFseq_vox( mri_smoothed, 0, 0, 0, 0 ), sizeof(float)*m_nDataSize);
    else
      memcpy(m_fData, &MRIFseq_vox( mri_smoothed, 0, 0, 0, 0 ), sizeof(float)*m_nDataSize);
    MRIfree(&mri_smoothed);
    MRIfree(&mri);
  }
  else
  {
    cerr << "Can not allocate mri\n";
    MRIfree(&mri);
  }
}

void SurfaceOverlay::SetActiveFrame(int nFrame)
{
  if (nFrame >= m_nNumOfFrames)
    nFrame = 0;
  m_nActiveFrame = nFrame;
  m_fData = m_fDataRaw + m_nActiveFrame*m_nDataSize;
  m_dMaxValue = m_dMinValue = m_fData[0];
  for ( int i = 0; i < m_nDataSize; i++ )
  {
    if ( m_dMaxValue < m_fData[i] )
    {
      m_dMaxValue = m_fData[i];
    }
    else if ( m_dMinValue > m_fData[i] )
    {
      m_dMinValue = m_fData[i];
    }
  }
  memcpy(m_fDataOriginal, m_fData, sizeof(float)*m_nDataSize);
}
