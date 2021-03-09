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


#include <assert.h>
#include "SurfaceOverlayProperty.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "FSSurface.h"
#include "SurfaceOverlay.h"
#include "SurfaceLabel.h"
#include <QDebug>
#include <QJsonDocument>
#include <QFile>

SurfaceOverlayProperty::SurfaceOverlayProperty ( SurfaceOverlay* overlay) :
  QObject( ),
  m_dOpacity( 1 ),
  m_bColorInverse( false ),
  m_bColorTruncate( false ),
  m_bClearLower(true),
  m_bClearHigher(false),
  m_bSmooth(false),
  m_nSmoothSteps(1),
  m_overlay( overlay ),
  m_bUsePercentile(false),
  m_dOffset(0),
  m_mask(NULL),
  m_maskData(NULL),
  m_bInverseMask(false),
  m_bIgnoreZeros(false),
  m_nColorScale(CS_Heat),
  m_nColorMethod(CM_LinearOpaque)
{
  m_lut = vtkRGBAColorTransferFunction::New();
}

SurfaceOverlayProperty::~SurfaceOverlayProperty ()
{
  m_lut->Delete();
  if (m_maskData)
    delete[] m_maskData;
  m_maskData = NULL;
}

void SurfaceOverlayProperty::Copy(SurfaceOverlayProperty *p)
{
  m_dOpacity = p->m_dOpacity;
  m_bColorInverse = p->m_bColorInverse;
  m_bColorTruncate = p->m_bColorTruncate;
  m_bClearLower = p->m_bClearLower;
  m_bClearHigher = p->m_bClearHigher;
  m_bSmooth = p->m_bSmooth;
  m_nSmoothSteps = p->m_nSmoothSteps;
  m_bUsePercentile = p->m_bUsePercentile;
  m_bIgnoreZeros = p->m_bIgnoreZeros;
  m_dOffset = p->m_dOffset;
  m_nColorScale = p->m_nColorScale;
  m_nColorMethod = p->m_nColorMethod;
  if (m_bUsePercentile)
  {
    m_dMinPoint = m_overlay->PercentileToPosition(p->m_overlay->PositionToPercentile(p->m_dMinPoint));
    m_dMidPoint = m_overlay->PercentileToPosition(p->m_overlay->PositionToPercentile(p->m_dMidPoint));
    m_dMaxPoint = m_overlay->PercentileToPosition(p->m_overlay->PositionToPercentile(p->m_dMaxPoint));
  }
  else
  {
    m_dMinPoint = p->m_dMinPoint;
    m_dMidPoint = p->m_dMidPoint;
    m_dMaxPoint = p->m_dMaxPoint;
  }
  m_customScale = p->m_customScale;
  m_dMinStop = p->m_dMinStop;
  m_dMaxStop = p->m_dMaxStop;
  for (int i = 0; i < 3; i++)
  {
    m_colorMin[i] = p->m_colorMin[i];
    m_colorMid[i] = p->m_colorMid[i];
    m_colorMax[i] = p->m_colorMax[i];
  }
  SetColorScale(m_nColorScale);
  m_bInverseMask = p->m_bInverseMask;
  m_mask = p->m_mask;
  if (p->m_maskData)
  {
    if (!m_maskData)
      m_maskData = new unsigned char[m_overlay->GetDataSize()];
    memcpy(m_maskData, p->m_maskData, m_overlay->GetDataSize());
  }
}

void SurfaceOverlayProperty::Reset()
{
  if ( m_overlay )
  {
    m_dMinPoint = m_overlay->PercentileToPosition(50);
    m_dMaxPoint = m_overlay->PercentileToPosition(99);
    if (m_dMinPoint < 0 && m_dMaxPoint > fabs(m_dMinPoint))
      m_dMinPoint = fabs(m_dMinPoint);
    m_dMidPoint = ( m_dMinPoint + m_dMaxPoint ) / 2;
    m_dOffset = 0;
    m_customScale.clear();
    m_customScale << QGradientStop(m_dMinPoint, Qt::red);
    m_dMinStop = m_dMinPoint;
    m_dMaxStop = m_dMaxPoint;
    SetColorScale( m_nColorScale );
    SetColorMethod( m_nColorMethod );
  }
}

double SurfaceOverlayProperty::GetOpacity() const
{
  return m_dOpacity;
}

void SurfaceOverlayProperty::SetOpacity( double opacity )
{
  if ( m_dOpacity != opacity )
  {
    m_dOpacity = opacity;
  }
}

int SurfaceOverlayProperty::GetColorScale() const
{
  return m_nColorScale;
}

void SurfaceOverlayProperty::SetColorScale( int nScale )
{
  m_nColorScale = nScale;
  switch ( nScale )
  {
  case CS_GreenRed:
    m_colorMin[0] = 0;
    m_colorMin[1] = 255;
    m_colorMin[2] = 0;
    m_colorMax[0] = 255;
    m_colorMax[1] = 0;
    m_colorMax[2] = 0;
    break;
  case CS_Heat:
    m_colorMin[0] = 255;
    m_colorMin[1] = 0;
    m_colorMin[2] = 0;
    m_colorMax[0] = 255;
    m_colorMax[1] = 255;
    m_colorMax[2] = 0;
    break;
  case CS_BlueRed:
    m_colorMin[0] = 0;
    m_colorMin[1] = 0;
    m_colorMin[2] = 255;
    m_colorMax[0] = 255;
    m_colorMax[1] = 0;
    m_colorMax[2] = 0;
    break;
  }

  m_lut->RemoveAllPoints();
  if ( nScale <= CS_BlueRed )
  {
    if ( !m_bColorTruncate || m_bColorInverse )
    {
      if ( m_bColorInverse )
      {
        m_lut->AddRGBAPoint( -m_dMaxPoint + m_dOffset, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
        m_lut->AddRGBAPoint( -m_dMidPoint + m_dOffset, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( -m_dMinPoint + m_dOffset, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( -m_dMinPoint + m_dOffset, 0.4, 0.4, 0.4, 1 );
        }
      }
      else
      {
        m_lut->AddRGBAPoint( -m_dMaxPoint + m_dOffset, 1-m_colorMax[0]/255.0, m_colorMax[1]/255.0, 1-m_colorMax[2]/255.0, 1 );
        m_lut->AddRGBAPoint( -m_dMidPoint + m_dOffset, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( -m_dMinPoint + m_dOffset, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( -m_dMinPoint + m_dOffset, 0.4, 0.4, 0.4, 1 );
        }
      }
      m_lut->AddRGBAPoint( -m_dMinPoint+ m_dOffset + 0.0000000001, 0.4, 0.4, 0.4, 1 );
    }

    m_lut->AddRGBAPoint( 0, 0.4, 0.4, 0.4, 1 );

    if ( ! (m_bColorTruncate && m_bColorInverse) )
    {
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset-0.0000000001, 0.4, 0.4, 0.4, 1 );
      if ( m_bColorInverse )
      {
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, 0.4, 0.4, 0.4, 1 );
        }
        m_lut->AddRGBAPoint( m_dMidPoint + m_dOffset, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        m_lut->AddRGBAPoint( m_dMaxPoint + m_dOffset, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
      }
      else
      {
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, 0.4, 0.4, 0.4, 1 );
        }
        m_lut->AddRGBAPoint( m_dMidPoint + m_dOffset, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        m_lut->AddRGBAPoint( m_dMaxPoint + m_dOffset, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
      }
    }
  }
  else if ( nScale == CS_ColorWheel)
  {
    if ( !m_bColorInverse )
    {
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, 1, 0, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.25, 1, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.5, 0, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.75, 0, 1, 1, 1);
      m_lut->AddRGBAPoint( m_dMaxPoint + m_dOffset, 0, 0, 1, 1);
    }
    else
    {
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset, 0, 0, 1, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.25, 0, 1, 1, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.5, 0, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + m_dOffset + (m_dMaxPoint-m_dMinPoint)*.75, 1, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMaxPoint + m_dOffset, 1, 0, 0, 1);
    }
  }
  else if ( nScale == CS_Custom)
  {
    for (int i = 0; i < m_customScale.size(); i++)
    {
      QColor c = m_customScale.at(i).second;
      m_lut->AddRGBAPoint(m_customScale.at(i).first + m_dOffset, c.redF(), c.greenF(), c.blueF(), 1);
    }
  }

  m_lut->Build();
}

void SurfaceOverlayProperty::SetCustomColorScale(QGradientStops stops)
{
  m_customScale = stops;
  m_dMinStop = stops[0].first;
  m_dMaxStop = m_dMinStop;
  for (int i = 0; i < stops.size(); i++)
  {
    if ( m_dMinStop > stops[i].first )
    {
      m_dMinStop = stops[i].first;
    }
    if (m_dMaxStop < stops[i].first )
    {
      m_dMaxStop = stops[i].first;
    }
  }
  SetColorScale(m_nColorScale);
}

int SurfaceOverlayProperty::GetColorMethod()
{
  return m_nColorMethod;
}

void SurfaceOverlayProperty::SetColorMethod( int nTh )
{
  m_nColorMethod = nTh;
}

void SurfaceOverlayProperty::SetMinPoint( double dValue )
{
  if ( dValue != m_dMinPoint )
  {
    m_dMinPoint = dValue;
    SetColorScale(m_nColorScale);
  }
}

double SurfaceOverlayProperty::GetMinPoint()
{
  return m_dMinPoint;
}


void SurfaceOverlayProperty::SetMidPoint( double dValue )
{
  if ( dValue != m_dMidPoint )
  {
    m_dMidPoint = dValue;
    SetColorScale(m_nColorScale);
  }
}

double SurfaceOverlayProperty::GetMidPoint()
{
  return m_dMidPoint;
}

void SurfaceOverlayProperty::SetMaxPoint( double dValue )
{
  if ( dValue != m_dMaxPoint )
  {
    m_dMaxPoint = dValue;
    SetColorScale(m_nColorScale);
  }
}

double SurfaceOverlayProperty::GetMaxPoint()
{
  return m_dMaxPoint;
}


void SurfaceOverlayProperty::SetOffset(double dOffset)
{
  if (dOffset != m_dOffset)
  {
    m_dOffset = dOffset;
    SetColorScale(m_nColorScale);
  }
}

double SurfaceOverlayProperty::GetOffset()
{
  return m_dOffset;
}

bool SurfaceOverlayProperty::GetColorInverse()
{
  return m_bColorInverse;
}

void SurfaceOverlayProperty::SetColorInverse( bool bInverse )
{
  m_bColorInverse = bInverse;
  SetColorScale( m_nColorScale );
}

bool SurfaceOverlayProperty::GetColorTruncate()
{
  return m_bColorTruncate;
}

void SurfaceOverlayProperty::SetColorTruncate( bool bTruncate )
{
  m_bColorTruncate = bTruncate;
  SetColorScale( m_nColorScale );
}

void SurfaceOverlayProperty::MapOverlayColor( float* data, unsigned char* colordata, int nPoints )
{
  if ( m_nColorScale <= CS_BlueRed )
  {
    MapOverlayColorSymmetric(data, colordata, nPoints);
  }
  else
  {
    MapOverlayColorFullScale(data, colordata, nPoints);
  }
}

void SurfaceOverlayProperty::MapOverlayColorSymmetric( float* data, unsigned char* colordata, int nPoints )
{
  double c[3];
  double dMidPoint = m_dMidPoint;
  if (m_nColorMethod != CM_Piecewise)
    dMidPoint = (m_dMaxPoint + m_dMinPoint)/2.0;
  if ( m_nColorMethod== CM_LinearOpaque )
  {
    for ( int i = 0; i < nPoints; i++ )
    {
      // map positive values
      bool bInMask = true;
      if (nPoints == m_overlay->GetDataSize())
        bInMask = (!m_mask || ( (!m_bInverseMask && m_maskData[i]) || (m_bInverseMask && !m_maskData[i]) ));
      if ( data[i] >= m_dMinPoint + m_dOffset && !( m_bColorInverse && m_bColorTruncate ) && bInMask )
      {
        double r = 0;
        if (m_dMaxPoint != dMidPoint)
        {
          r = ( m_dMaxPoint + m_dOffset - data[i] ) / ( m_dMaxPoint - dMidPoint );
        }
        if ( r < 0 )
        {
          r = 0;
        }
        else if ( r > 1 )
        {
          r = 1;
        }

        if ( !m_bColorInverse )
        {
          c[0] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
          c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
          c[2] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
        }
        else
        {
          c[0] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
          c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
          c[2] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
        }
        colordata[i*4]    = ( int )( colordata[i*4]   * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1]  = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2]  = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
      // map negative values
      else if ( data[i] <= -m_dMinPoint + m_dOffset && !( m_bColorTruncate && !m_bColorInverse ) && bInMask )
      {
        double r = 0;
        if (m_dMaxPoint != dMidPoint)
        {
          r = ( data[i] + m_dMaxPoint + m_dOffset ) / ( m_dMaxPoint - dMidPoint );
        }
        if ( r < 0 )
        {
          r = 0;
        }
        else if ( r > 1 )
        {
          r = 1;
        }

        if ( !m_bColorInverse )
        {
          c[0] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
          c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
          c[2] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
        }
        else
        {
          c[0] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
          c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
          c[2] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
        }
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
    }
  }
  else if ( m_nColorMethod == CM_Piecewise || m_nColorMethod == CM_Linear )
  {
    for ( int i = 0; i < nPoints; i++ )
    {
      // map positive values
      bool bInMask = true;
      if (nPoints == m_overlay->GetDataSize())
        bInMask = (!m_mask || ( (!m_bInverseMask && m_maskData[i]) || (m_bInverseMask && !m_maskData[i]) ));
      if ( data[i] >= m_dMinPoint + m_dOffset && !( m_bColorInverse && m_bColorTruncate ) && bInMask )
      {
        if ( data[i] < dMidPoint + m_dOffset )
        {
          double r = 0;
          if (dMidPoint != m_dMinPoint)
          {
            r = ( dMidPoint + m_dOffset - data[i] ) / ( dMidPoint - m_dMinPoint );
          }
          if ( r < 0 )
          {
            r = 0;
          }
          else if ( r > 1 )
          {
            r = 1;
          }
          if ( !m_bColorInverse )
          {
            c[0] = colordata[i*4] * r + m_colorMin[0] * ( 1 - r );
            c[1] = colordata[i*4+1] * r + m_colorMin[1] * ( 1 - r );
            c[2] = colordata[i*4+2] * r + m_colorMin[2] * ( 1 - r );
          }
          else
          {
            c[0] = colordata[i*4] * r + m_colorMin[2] * ( 1 - r );
            c[1] = colordata[i*4+1] * r + m_colorMin[1] * ( 1 - r );
            c[2] = colordata[i*4+2] * r + m_colorMin[0] * ( 1 - r );
          }
        }
        else if ( data[i] <= m_dMaxPoint + m_dOffset )
        {
          double r = 0;
          if (m_dMaxPoint != dMidPoint)
          {
            r = ( m_dMaxPoint + m_dOffset - data[i] ) / ( m_dMaxPoint - dMidPoint );
          }
          if ( r < 0 )
          {
            r = 0;
          }
          else if ( r > 1 )
          {
            r = 1;
          }
          if ( !m_bColorInverse )
          {
            c[0] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
            c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
            c[2] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
          }
          else
          {
            c[0] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
            c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
            c[2] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
          }
        }
        else
        {
          if ( !m_bColorInverse )
          {
            c[0] = m_colorMax[0];
            c[1] = m_colorMax[1];
            c[2] = m_colorMax[2];
          }
          else
          {
            c[0] = m_colorMax[2];
            c[1] = m_colorMax[1];
            c[2] = m_colorMax[0];
          }
        }
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
      // map negative value
      else if ( data[i] <= -m_dMinPoint + m_dOffset && !( m_bColorTruncate && !m_bColorInverse ) && bInMask )
      {
        if ( data[i] >= -dMidPoint + m_dOffset )
        {
          double r = 0;
          if (m_dMinPoint != dMidPoint)
          {
            r = ( dMidPoint + m_dOffset + data[i] ) / ( dMidPoint - m_dMinPoint );
          }
          if ( r < 0 )
          {
            r = 0;
          }
          else if ( r > 1 )
          {
            r = 1;
          }
          if ( !m_bColorInverse )
          {
            c[0] = colordata[i*4] * r + m_colorMin[2] * ( 1 - r );
            c[1] = colordata[i*4+1] * r + m_colorMin[1] * ( 1 - r );
            c[2] = colordata[i*4+2] * r + m_colorMin[0] * ( 1 - r );
          }
          else
          {
            c[0] = colordata[i*4] * r + m_colorMin[0] * ( 1 - r );
            c[1] = colordata[i*4+1] * r + m_colorMin[1] * ( 1 - r );
            c[2] = colordata[i*4+2] * r + m_colorMin[2] * ( 1 - r );
          }
        }
        else if ( data[i] >= -m_dMaxPoint + m_dOffset )
        {
          double r = 0;
          if (m_dMaxPoint != dMidPoint)
          {
            r = ( m_dMaxPoint + m_dOffset + data[i] ) / ( m_dMaxPoint - dMidPoint );
          }
          if ( r < 0 )
          {
            r = 0;
          }
          else if ( r > 1 )
          {
            r = 1;
          }
          if ( !m_bColorInverse )
          {
            c[0] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
            c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
            c[2] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
          }
          else
          {
            c[0] = m_colorMin[0] * r + m_colorMax[0] * ( 1 - r );
            c[1] = m_colorMin[1] * r + m_colorMax[1] * ( 1 - r );
            c[2] = m_colorMin[2] * r + m_colorMax[2] * ( 1 - r );
          }
        }
        else
        {
          if ( !m_bColorInverse )
          {
            c[0] = m_colorMax[2];
            c[1] = m_colorMax[1];
            c[2] = m_colorMax[0];
          }
          else
          {
            c[0] = m_colorMax[0];
            c[1] = m_colorMax[1];
            c[2] = m_colorMax[2];
          }
        }
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
    }
  }
}

void SurfaceOverlayProperty::MapOverlayColorFullScale( float* data, unsigned char* colordata, int nPoints )
{
  if (!m_overlay)
  {
    LABEL* label = NULL;
    if (m_mask)
      label = m_mask->GetLabelData();
    return;
  }
  double c[4];
  double dThLow = m_dMinPoint + m_dOffset;
  double dThHigh = m_overlay->m_dMaxValue+1e10;
  if ( m_nColorScale == CS_Custom)
  {
    if ( m_bClearLower )
    {
      dThLow = m_dMinStop + m_dOffset;
    }
    else
    {
      dThLow = m_overlay->m_dMinValue-1e10;
    }
    if (m_bClearHigher)
    {
      dThHigh = m_dMaxStop + m_dOffset;
    }
  }
  for ( int i = 0; i < nPoints; i++ )
  {
    bool bInMask = true;
    if (nPoints == m_overlay->GetDataSize())
      bInMask = (!m_mask || ( (!m_bInverseMask && m_maskData[i]) || (m_bInverseMask && !m_maskData[i]) ));
    if (bInMask && data[i] >= dThLow && data[i] <= dThHigh)
    {
      m_lut->GetColor( data[i], c );
      colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0]*255 * m_dOpacity );
      colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1]*255 * m_dOpacity );
      colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2]*255 * m_dOpacity );
    }
  }
}

void SurfaceOverlayProperty::SetClearLower(bool bClear)
{
  m_bClearLower = bClear;
  SetColorScale(m_nColorScale);
}

void SurfaceOverlayProperty::SetClearHigher(bool bClear)
{
  m_bClearHigher = bClear;
  SetColorScale(m_nColorScale);
}

void SurfaceOverlayProperty::SetSmooth(bool bSmooth)
{
  if (m_bSmooth != bSmooth)
  {
    m_bSmooth = bSmooth;
    emit SmoothChanged();
  }
}

void SurfaceOverlayProperty::SetSmoothSteps(int n)
{
  if (n != m_nSmoothSteps)
  {
    m_nSmoothSteps = n;
    emit SmoothChanged();
  }
}

void SurfaceOverlayProperty::SetMask(SurfaceLabel *label)
{
  if (label != m_mask)
  {
    m_mask = label;
    if (!m_maskData)
      m_maskData = new unsigned char[m_overlay->GetDataSize()];
    memset(m_maskData, 0, m_overlay->GetDataSize());
    if (label)
    {
      LABEL* l = label->GetLabelData();
      for (int i = 0; i < l->n_points; i++)
        m_maskData[l->lv[i].vno] = 1;
      connect(label, SIGNAL(destroyed(QObject*)), SLOT(OnLabelMaskDestroyed(QObject*)), Qt::UniqueConnection);
    }
    emit MaskChanged();
  }
}

void SurfaceOverlayProperty::SetMaskInverse(bool b)
{
  if (b != m_bInverseMask)
  {
    m_bInverseMask = b;
    emit MaskChanged();
  }
}

void SurfaceOverlayProperty::OnLabelMaskDestroyed(QObject* label)
{
  if (label == m_mask)
  {
    SetMask(NULL);
    emit ColorMapChanged();
  }
}

bool SurfaceOverlayProperty::LoadCustomColorScale(const QString &filename)
{
  QFile file(filename);
  file.open(QIODevice::ReadOnly);
  QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
  QVariantList list = doc.toVariant().toList();
  file.close();
  if (list.isEmpty())
  {
    qDebug() << "Unable to load color scale from " << filename;
    return false;
  }
  else
  {
    QGradientStops stops;
    foreach (QVariant v, list)
    {
      QVariantMap map = v.toMap();
      stops << QGradientStop(map["val"].toDouble(), QColor(map["r"].toInt(), map["g"].toInt(), map["b"].toInt()));
    }
    SetColorScale(CS_Custom);
    SetCustomColorScale(stops);
    return true;
  }
}

bool SurfaceOverlayProperty::SaveCustomColorScale(const QString &filename)
{
  QFile file(filename);
  if (file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    QVariantList list;
    foreach (QGradientStop stop, m_customScale)
    {
      QVariantMap map;
      map["val"] = stop.first;
      map["r"] = stop.second.red();
      map["g"] = stop.second.green();
      map["b"] = stop.second.blue();
      list << map;
    }
    file.write(QJsonDocument::fromVariant(list).toJson());
    return true;
  }
  else
  {
    qDebug() << "Unable to save color scale to " << filename;
    return false;
  }
}
