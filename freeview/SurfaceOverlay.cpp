/**
 * @brief Implementation for surface layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "LayerMRI.h"
#include "MyUtils.h"
#include <QDateTime>
#include "utils.h"

SurfaceOverlay::SurfaceOverlay ( LayerSurface* surf ) :
  QObject(),
  m_fData( NULL ),
  m_fDataRaw( NULL ),
  m_fDataUnsmoothed(NULL),
  m_dMaxValue(0),
  m_dMinValue(0),
  m_dNonZeroMinValue(0),
  m_dRawMaxValue(0),
  m_dRawMinValue(0),
  m_surface( surf ),
  m_bCorrelationData( false ),
  m_mriCorrelation(NULL),
  m_overlayPaired(NULL),
  m_nActiveFrame(0),
  m_nNumOfFrames(1),
  m_bComputeCorrelation(false),
  m_volumeCorrelationSource(NULL),
  m_fCorrelationSourceData(NULL),
  m_fCorrelationDataBuffer(NULL),
  m_ctab(NULL)
{
  InitializeData();

  m_nID = QDateTime::currentMSecsSinceEpoch();

  m_property =  new SurfaceOverlayProperty( this );
  connect( m_property, SIGNAL(ColorMapChanged()), surf, SLOT(UpdateOverlay()), Qt::UniqueConnection);
  connect( m_property, SIGNAL(SmoothChanged()), this, SLOT(UpdateSmooth()), Qt::UniqueConnection);
}

SurfaceOverlay::~SurfaceOverlay ()
{
  if ( m_fDataRaw )
    delete[] m_fDataRaw;

  if (m_fData)
    delete[] m_fData;

  if (m_fDataUnsmoothed)
    delete[] m_fDataUnsmoothed;

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

  if (m_fCorrelationSourceData)
    delete[] m_fCorrelationSourceData;

  if (m_fCorrelationDataBuffer)
    delete[] m_fCorrelationDataBuffer;

  if (m_ctab)
    ::CTABfree(&m_ctab);
}

void SurfaceOverlay::InitializeData()
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    if ( m_fDataRaw )
      delete[] m_fDataRaw;

    m_nDataSize = mris->nvertices;
    m_fDataRaw = new float[ m_nDataSize ];
    if ( !m_fDataRaw )
    {
      return;
    }

    if ( m_fData )
      delete[] m_fData;

    m_fData = new float[m_nDataSize];
    if ( !m_fData )
    {
      return;
    }

    m_fDataUnsmoothed = new float[m_nDataSize];
    if ( !m_fDataUnsmoothed )
    {
      return;
    }

    m_dMaxValue = m_dMinValue = mris->vertices[0].val;
    m_dNonZeroMinValue = 1e10;
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
      if (m_fData[vno] > 0 && m_dNonZeroMinValue > m_fData[vno])
        m_dNonZeroMinValue = m_fData[vno];
    }
    m_dRawMaxValue = m_dMaxValue;
    m_dRawMinValue = m_dMinValue;
    m_dDisplayRange[0] = m_dMinValue;
    m_dDisplayRange[1] = m_dMaxValue;
    memcpy(m_fDataRaw, m_fData, sizeof(float)*m_nDataSize);
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
    m_fDataRaw = data_buffer_in;
    if ( !m_fDataRaw )
      return;

    m_dMaxValue = m_dMinValue = m_fDataRaw[0];
    m_dNonZeroMinValue = 1e10;
    for (int i = 0; i < m_nDataSize*m_nNumOfFrames; i++)
    {
      if ( m_dMaxValue < m_fDataRaw[i])
        m_dMaxValue = m_fDataRaw[i];
      else if (m_dMinValue > m_fDataRaw[i])
        m_dMinValue = m_fDataRaw[i];
      if ( m_fDataRaw[i] > 0 && m_dNonZeroMinValue > m_fDataRaw[i])
        m_dNonZeroMinValue = m_fDataRaw[i];
    }
    m_dRawMaxValue = m_dMaxValue;
    m_dRawMinValue = m_dMinValue;

    m_dDisplayRange[0] = m_dMinValue;
    m_dDisplayRange[1] = m_dMaxValue;

    if ( m_fData )
      delete[] m_fData;

    m_fData = new float[ m_nDataSize ];
    if ( !m_fData )
      return;

    m_fDataUnsmoothed = new float[m_nDataSize];
    if (!m_fDataUnsmoothed)
      return;

    SetActiveFrame(0);
    GetProperty()->Reset();
    m_fCorrelationSourceData = new float[nframes];
    memset(m_fCorrelationSourceData, 0, sizeof(float)*nframes);
    m_fCorrelationDataBuffer = new float[nframes];
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
  MRI* mri = ::MRIreadHeader( filename.toLatin1().data(), -1 );
  if ( mri == NULL )
  {
    cerr << "MRIread failed: unable to read from " << qPrintable(filename) << "\n";
    return false;
  }
  if ( (((qlonglong)mri->width)*mri->height*mri->nframes)%(m_nDataSize*m_nDataSize) != 0 )
  {
    cerr << "Correlation data does not match with surface\n";
    MRIfree( &mri );
    return false;
  }
  MRIfree( &mri );
  ::SetProgressCallback(ProgressCallback, 0, 100);
  try {
    mri = ::MRIread( filename.toLatin1().data() );      // long process
  }
  catch (int ret) {
    return false;
  }

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
  bool bSingleHemiData = ( ((qlonglong)m_mriCorrelation->width)*m_mriCorrelation->height*m_mriCorrelation->nframes == m_nDataSize*m_nDataSize );
  if ( nHemisphere == -1)
  {
    nHemisphere = m_surface->GetHemisphere();
  }
  int nVertexOffset = nHemisphere * m_nDataSize;
  int nDataOffset = m_surface->GetHemisphere() * m_nDataSize;
  if (bSingleHemiData)
  {
    nVertexOffset = 0;
    nDataOffset = 0;
  }
  double old_range = m_dMaxValue - m_dMinValue;
  if (m_mriCorrelation->height > 1)
    m_dMaxValue = m_dMinValue =
        MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, nDataOffset, 0, 0 );
  else
    m_dMaxValue = m_dMinValue =
        MRIFseq_vox( m_mriCorrelation, nVertex + nVertexOffset, 0, 0, nDataOffset );
  m_dNonZeroMinValue = 1e10;
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
    if (m_fData[i] > 0 && m_dNonZeroMinValue > m_fData[i])
      m_dNonZeroMinValue = m_fData[i];
  }
  memcpy(m_fDataRaw, m_fData, sizeof(float)*m_nDataSize);
  memcpy(m_fDataUnsmoothed, m_fData, sizeof(float)*m_nDataSize);
  if (GetProperty()->GetSmooth())
  {
    SmoothData();
  }

  m_bCorrelationDataReady = true;
  if (old_range <= 0)
  {
    m_property->Reset();
  }

  if (m_overlayPaired && nHemisphere == m_surface->GetHemisphere() && !bSingleHemiData)
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
    memcpy(m_fData, m_fDataUnsmoothed, sizeof(float)*m_nDataSize);
  }
  m_surface->UpdateOverlay(true);
  emit DataUpdated();

  if (trigger_paired && m_overlayPaired)
    m_overlayPaired->UpdateSmooth(false);
}

void SurfaceOverlay::SmoothData(int nSteps_in, float *data_out)
{
  MRI* mri = NULL;
  try {
    mri = MRIallocSequence( m_nDataSize,
                            1,
                            1,
                            MRI_FLOAT, 1);
  } catch (int ret) {
    return;
  }

  if (!mri)
  {
    cerr << "Can not allocate mri\n";
    return;
  }
  memcpy(&MRIFseq_vox( mri, 0, 0, 0, 0 ), m_fDataUnsmoothed, sizeof(float)*m_nDataSize);
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
  memcpy(m_fData, m_fDataRaw + m_nActiveFrame*m_nDataSize, sizeof(float)*m_nDataSize);
  memcpy(m_fDataUnsmoothed, m_fData, sizeof(float)*m_nDataSize);
  m_dMaxValue = m_dMinValue = m_fData[0];
  m_dNonZeroMinValue = 1e10;
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
    if (m_fData[i] > 0 && m_dNonZeroMinValue > m_fData[i])
      m_dNonZeroMinValue = m_fData[i];
  }
  if (GetProperty()->GetSmooth())
  {
    SmoothData();
  }
}

void SurfaceOverlay::SetComputeCorrelation(bool flag)
{
  m_bComputeCorrelation = flag;
  if (flag)
  {
    UpdateCorrelationCoefficient();
  }
  else
  {
    SetActiveFrame(m_nActiveFrame);
  }
}

void SurfaceOverlay::UpdateCorrelationCoefficient(double* pos_in)
{
  if (m_bComputeCorrelation)
  {
    if (m_volumeCorrelationSource && m_volumeCorrelationSource->GetNumberOfFrames() == m_nNumOfFrames)
    {
      double pos[3];
      int n[3];
      m_volumeCorrelationSource->GetSlicePosition(pos);
      m_volumeCorrelationSource->TargetToRAS(pos, pos);
      m_volumeCorrelationSource->RASToOriginalIndex(pos, n);
      m_volumeCorrelationSource->GetVoxelValueByOriginalIndexAllFrames(n[0], n[1], n[2], m_fCorrelationSourceData);
    }
    else
    {
      int nVertex = m_surface->GetCurrentVertex();
      if (!m_surface->IsInflated() && pos_in)
        nVertex = m_surface->GetVertexIndexAtTarget(pos_in, NULL);

      if (nVertex >= 0)
      {
        this->GetDataAtVertex(nVertex, m_fCorrelationSourceData);
      }
    }
    for (int i = 0; i < m_nDataSize; i++)
    {
      for (int j = 0; j < m_nNumOfFrames; j++)
      {
        m_fCorrelationDataBuffer[j] = m_fDataRaw[i+j*m_nDataSize];
      }
      m_fData[i] = MyUtils::CalculateCorrelationCoefficient(m_fCorrelationSourceData, m_fCorrelationDataBuffer, m_nNumOfFrames);

      if (i == 0)
      {
        m_dMaxValue = m_dMinValue = m_fData[0];
        if (m_fData[0] > 0)
          m_dNonZeroMinValue = m_fData[0];
        else
          m_dNonZeroMinValue = 1e10;
      }
      else if ( m_dMaxValue < m_fData[i] )
      {
        m_dMaxValue = m_fData[i];
      }
      else if ( m_dMinValue > m_fData[i] )
      {
        m_dMinValue = m_fData[i];
      }
      if (m_fData[i] > 0 && m_dNonZeroMinValue > m_fData[i])
        m_dNonZeroMinValue = m_fData[i];
    }
    memcpy(m_fDataRaw, m_fData, sizeof(float)*m_nDataSize);
    memcpy(m_fDataUnsmoothed, m_fData, sizeof(float)*m_nDataSize);
    if (GetProperty()->GetSmooth())
      SmoothData();
    m_surface->UpdateOverlay(true);
    emit DataUpdated();
  }
}

void SurfaceOverlay::SetCorrelationSourceVolume(LayerMRI *vol)
{
  if (m_volumeCorrelationSource)
    m_volumeCorrelationSource->disconnect(this, 0);
  m_volumeCorrelationSource = vol;
  if (vol)
    connect(vol, SIGNAL(destroyed(QObject*)), this, SLOT(OnCorrelationSourceDeleted(QObject*)));
  UpdateCorrelationCoefficient();
}

void SurfaceOverlay::OnCorrelationSourceDeleted(QObject *obj)
{
  if (m_volumeCorrelationSource == obj)
  {
    m_volumeCorrelationSource = NULL;
    UpdateCorrelationCoefficient();
  }
}

void SurfaceOverlay::GetRange( double* range )
{
  if (m_bComputeCorrelation)
  {
    range[0] = -1;
    range[1] = 1;
  }
  else
  {
    range[0] = m_dMinValue;
    range[1] = m_dMaxValue;
  }
}

void SurfaceOverlay::GetNonZeroRange(double *range)
{
  if (m_bComputeCorrelation)
  {
    range[0] = -1;
    range[1] = 1;
  }
  else
  {
    range[0] = m_dNonZeroMinValue;
    range[1] = m_dMaxValue;
  }
}

void SurfaceOverlay::GetRawRange( double* range )
{
  range[0] = m_dRawMinValue;
  range[1] = m_dRawMaxValue;
}

bool SurfaceOverlay::GetDataAtVertex(int nVertex, float *output)
{
  if (nVertex < 0)
    return false;

  for (int i = 0; i < m_nNumOfFrames; i++)
  {
    output[i] = m_fDataRaw[nVertex+i*m_nDataSize];
  }
  return true;
}

double SurfaceOverlay::PercentileToPosition(double dPercentile)
{
  return PercentileToPosition(dPercentile, GetProperty()->GetIgnoreZeros());
}

double SurfaceOverlay::PercentileToPosition(double percentile_in, bool bIgnoreZeros)
{
  double percentile = percentile_in/100;
  double range[2];
  if (bIgnoreZeros)
    GetNonZeroRange(range);
  else
    GetRange(range);
  int m_nNumberOfBins = 100;
  double m_dBinWidth = ( range[1] - range[0] ) / m_nNumberOfBins;
  int* m_nOutputData = new int[m_nNumberOfBins];
  if ( !m_nOutputData )
  {
    qCritical() << "Can not allocate memory.";
    return 0;
  }

  // calculate histogram data
  memset( m_nOutputData, 0, m_nNumberOfBins * sizeof( int ) );
  for ( long i = 0; i < m_nDataSize; i++ )
  {
    if (!bIgnoreZeros || m_fData[i] != 0)
    {
      int n = (int)( ( m_fData[i] - range[0] ) / m_dBinWidth );
      if ( n >= 0 && n < m_nNumberOfBins )
      {
        m_nOutputData[n] ++;
      }
    }
  }

  double m_dOutputTotalArea = 0;
  for (int i = 0; i < m_nNumberOfBins; i++)
  {
    m_dOutputTotalArea += m_nOutputData[i];
  }

  double dArea = 0;
  double dPos = range[0];
  int n = 0;
  while (dArea/m_dOutputTotalArea < percentile && n < m_nNumberOfBins)
  {
    dArea += m_nOutputData[n];
    dPos += m_dBinWidth;
    n++;
  }

  if (dArea > percentile*m_dOutputTotalArea && n > 0)
  {
    dPos -= (dArea-percentile*m_dOutputTotalArea)*m_dBinWidth/m_nOutputData[n-1];
  }

  delete[] m_nOutputData;
  return dPos;
}

double SurfaceOverlay::PositionToPercentile(double dPos)
{
  return PositionToPercentile(dPos, GetProperty()->GetIgnoreZeros());
}

double SurfaceOverlay::PositionToPercentile(double pos, bool bIgnoreZeros)
{
  double range[2];
  if (bIgnoreZeros)
    GetNonZeroRange(range);
  else
    GetRange(range);
  int m_nNumberOfBins = 100;
  double m_dBinWidth = ( range[1] - range[0] ) / m_nNumberOfBins;
  int* m_nOutputData = new int[m_nNumberOfBins];
  if ( !m_nOutputData )
  {
    qCritical() << "Can not allocate memory.";
    return 0;
  }

  // calculate histogram data
  memset( m_nOutputData, 0, m_nNumberOfBins * sizeof( int ) );
  for ( long i = 0; i < m_nDataSize; i++ )
  {
    if (!bIgnoreZeros || m_fData[i] != 0)
    {
      int n = (int)( ( m_fData[i] - range[0] ) / m_dBinWidth );
      if ( n >= 0 && n < m_nNumberOfBins )
      {
        m_nOutputData[n] ++;
      }
    }
  }

  double m_dOutputTotalArea = 0;
  for (int i = 0; i < m_nNumberOfBins; i++)
  {
    m_dOutputTotalArea += m_nOutputData[i];
  }

  double dArea = 0;
  double dPos = range[0];
  int n = 0;
  while (dPos < pos && n < m_nNumberOfBins)
  {
    dArea += m_nOutputData[n];
    dPos += m_dBinWidth;
    n++;
  }
  if (dPos > pos && n > 0)
  {
    dArea -= (dPos-pos)*m_nOutputData[n-1]/m_dBinWidth;
  }

  return 100*dArea/m_dOutputTotalArea;
}

void SurfaceOverlay::UpdateMaxHistCount(double* range, int nBins)
{
  int* CntData = new int[nBins];
  memset( CntData, 0, nBins * sizeof( int ) );
  double binWidth = ( range[1] - range[0] ) / nBins;
  for ( long i = 0; i < m_nDataSize; i++ )
  {
    int n = (int)( ( m_fData[i] - range[0] ) / binWidth );
    if ( n >= 0 && n < nBins )
    {
      CntData[n] ++;
    }
  }

  // find max and second max
  int nMaxCount = 0;
  for ( int i = 0; i < nBins; i++ )
  {
    if ( nMaxCount < CntData[i] )
    {
      nMaxCount  = CntData[i];
    }
  }

  setProperty("HistMaxCount", nMaxCount);
  setProperty("HistRange", range[0]);
  setProperty("HistBins", nBins);
}

void SurfaceOverlay::SetColorTable(COLOR_TABLE *ctab)
{
  if (ctab)
    m_ctab = ::CTABdeepCopy(ctab);
}
