/**
 * @file  SurfaceOverlayProperty.cxx
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
 *    $Date: 2012/06/26 17:06:07 $
 *    $Revision: 1.6 $
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


#include <assert.h>
#include "SurfaceOverlayProperty.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "FSSurface.h"
#include "SurfaceOverlay.h"
#include <QDebug>

SurfaceOverlayProperty::SurfaceOverlayProperty ( SurfaceOverlay* overlay) :
  QObject( ),
  m_dOpacity( 1 ),
  m_bColorInverse( false ),
  m_bColorTruncate( false ),
  m_bClearLower(true),
  m_bClearHigher(false),
  m_bSmooth(false),
  m_nSmoothSteps(1),
  m_overlay( overlay )
{
  m_lut = vtkRGBAColorTransferFunction::New();

  Reset();
  SetColorScale( CS_Heat );
  SetColorMethod( CM_LinearOpaque );
}

SurfaceOverlayProperty::~SurfaceOverlayProperty ()
{
  m_lut->Delete();
}

void SurfaceOverlayProperty::Reset()
{
  if ( m_overlay )
  {
    m_dMinPoint = fabs( m_overlay->m_dMinValue + m_overlay->m_dMaxValue ) / 2;
    m_dMaxPoint = m_overlay->m_dMaxValue;
    m_dMidPoint = ( m_dMinPoint + m_dMaxPoint ) / 2;
    m_customScale.clear();
    m_customScale << QGradientStop(m_dMinPoint, Qt::red) << QGradientStop(m_dMaxPoint, Qt::yellow);
    m_dMinStop = m_dMinPoint;
    m_dMaxStop = m_dMaxPoint;
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
        m_lut->AddRGBAPoint( -m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
        m_lut->AddRGBAPoint( -m_dMidPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( -m_dMinPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( -m_dMinPoint, 0.4, 0.4, 0.4, 1 );
        }
      }
      else
      {
        m_lut->AddRGBAPoint( -m_dMaxPoint, 1-m_colorMax[0]/255.0, m_colorMax[1]/255.0, 1-m_colorMax[2]/255.0, 1 );
        m_lut->AddRGBAPoint( -m_dMidPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( -m_dMinPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( -m_dMinPoint, 0.4, 0.4, 0.4, 1 );
        }
      }
      m_lut->AddRGBAPoint( -m_dMinPoint+0.0000000001, 0.4, 0.4, 0.4, 1 );
    }

    m_lut->AddRGBAPoint( 0, 0.4, 0.4, 0.4, 1 );

    if ( ! (m_bColorTruncate && m_bColorInverse) )
    {
      m_lut->AddRGBAPoint( m_dMinPoint-0.0000000001, 0.4, 0.4, 0.4, 1 );
      if ( m_bColorInverse )
      {
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( m_dMinPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( m_dMinPoint, 0.4, 0.4, 0.4, 1 );
        }
        m_lut->AddRGBAPoint( m_dMidPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
        m_lut->AddRGBAPoint( m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
      }
      else
      {
        if ( m_nColorMethod == CM_LinearOpaque )
        {
          m_lut->AddRGBAPoint( m_dMinPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        }
        else
        {
          m_lut->AddRGBAPoint( m_dMinPoint, 0.4, 0.4, 0.4, 1 );
        }
        m_lut->AddRGBAPoint( m_dMidPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
        m_lut->AddRGBAPoint( m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
      }
    }
  }
  else if ( nScale == CS_ColorWheel)
  {
    if ( !m_bColorInverse )
    {
      m_lut->AddRGBAPoint( m_dMinPoint, 1, 0, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.25, 1, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.5, 0, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.75, 0, 1, 1, 1);
      m_lut->AddRGBAPoint( m_dMaxPoint, 0, 0, 1, 1);
    }
    else
    {
      m_lut->AddRGBAPoint( m_dMinPoint, 0, 0, 1, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.25, 0, 1, 1, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.5, 0, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMinPoint + (m_dMaxPoint-m_dMinPoint)*.75, 1, 1, 0, 1);
      m_lut->AddRGBAPoint( m_dMaxPoint, 1, 0, 0, 1);
    }
  }
  else if ( nScale == CS_Custom)
  {
    for (int i = 0; i < m_customScale.size(); i++)
    {
      QColor c = m_customScale.at(i).second;
      m_lut->AddRGBAPoint(m_customScale.at(i).first, c.redF(), c.greenF(), c.blueF(), 1);
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
  }
}

double SurfaceOverlayProperty::GetMaxPoint()
{
  return m_dMaxPoint;
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

/*
void SurfaceOverlayProperty::MapOverlayColor( unsigned char* colordata, int nPoints )
{
  MapOverlayColor( m_overlay->GetData(), colordata, nPoints );
}
*/

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
      if ( data[i] >= m_dMinPoint && !( m_bColorInverse && m_bColorTruncate ) )
      {
        double r = 0;
        if (m_dMaxPoint != dMidPoint)
        {
          r = ( m_dMaxPoint - data[i] ) / ( m_dMaxPoint - dMidPoint );
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
      else if ( data[i] <= -m_dMinPoint && !( m_bColorTruncate && !m_bColorInverse ) )
      {
        double r = 0;
        if (m_dMaxPoint != dMidPoint)
        {
          r = ( data[i] + m_dMaxPoint ) / ( m_dMaxPoint - dMidPoint );
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
      if ( data[i] >= m_dMinPoint && !( m_bColorInverse && m_bColorTruncate ) )
      {
        if ( data[i] < dMidPoint )
        {
          double r = 0;
          if (dMidPoint != m_dMinPoint)
          {
            r = ( dMidPoint - data[i] ) / ( dMidPoint - m_dMinPoint );
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
        else if ( data[i] <= m_dMaxPoint )
        {
          double r = 0;
          if (m_dMaxPoint != dMidPoint)
          {
            r = ( m_dMaxPoint - data[i] ) / ( m_dMaxPoint - dMidPoint );
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
      else if ( data[i] <= -m_dMinPoint && !( m_bColorTruncate && !m_bColorInverse ) )
      {
        if ( data[i] >= -dMidPoint )
        {
          double r = 0;
          if (m_dMinPoint != dMidPoint)
          {
            r = ( dMidPoint + data[i] ) / ( dMidPoint - m_dMinPoint );
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
        else if ( data[i] >= -m_dMaxPoint )
        {
          double r = 0;
          if (m_dMaxPoint != dMidPoint)
          {
            r = ( m_dMaxPoint + data[i] ) / ( m_dMaxPoint - dMidPoint );
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
    return;
  }
  double c[4];
  double dThLow = m_dMinPoint;
  double dThHigh = m_overlay->m_dMaxValue+1e10;
  if ( m_nColorScale == CS_Custom)
  {
    if ( m_bClearLower )
    {
      dThLow = m_dMinStop;
    }
    else
    {
      dThLow = m_overlay->m_dMinValue-1e10;
    }
    if (m_bClearHigher)
    {
      dThHigh = m_dMaxStop;
    }
  }
  for ( int i = 0; i < nPoints; i++ )
  {
    if (data[i] >= dThLow && data[i] <= dThHigh)
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
