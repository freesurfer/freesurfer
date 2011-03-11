/**
 * @file  SurfaceOverlayProperties.cxx
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
#include "SurfaceOverlayProperties.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "FSSurface.h"
#include "SurfaceOverlay.h"
#include <wx/filename.h>

SurfaceOverlayProperties::SurfaceOverlayProperties ( SurfaceOverlay* overlay) :
    Broadcaster( "SurfaceOverlayProperties" ),
    m_dOpacity( 1 ),
    m_bColorInverse( false ),
    m_bColorTruncate( false ),
    m_overlay( overlay )
{ 
  m_lut = vtkRGBAColorTransferFunction::New();
  
  if ( overlay )
  {
    m_dMinPoint = fabs( overlay->m_dMinValue + overlay->m_dMaxValue ) / 2;
    m_dMaxPoint = overlay->m_dMaxValue;
    m_dMidPoint = ( m_dMinPoint + m_dMaxPoint ) / 2;
  }
  SetColorScale( CS_Heat );
  SetColorMethod( CM_LinearOpaque );
}

SurfaceOverlayProperties::~SurfaceOverlayProperties ()
{
  m_lut->Delete();
}

void SurfaceOverlayProperties::ColorMapChanged()
{
  // Notify the layers that use the color map stuff.
  this->SendBroadcast( "ColorMapChanged", this );
}

void SurfaceOverlayProperties::SetSurfaceOverlay( SurfaceOverlay* overlay )
{
  m_overlay = overlay;

  this->ColorMapChanged();
}

double SurfaceOverlayProperties::GetOpacity() const
{
  return m_dOpacity;
}

void SurfaceOverlayProperties::SetOpacity( double opacity )
{
  if ( m_dOpacity != opacity )
  {
    m_dOpacity = opacity;
  }
}

int SurfaceOverlayProperties::GetColorScale() const
{
  return m_nColorScale;
}

void SurfaceOverlayProperties::SetColorScale( int nScale )
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
  if ( !m_bColorTruncate || m_bColorInverse )
  {
    if ( m_bColorInverse )
    {
      m_lut->AddRGBAPoint( -m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
      m_lut->AddRGBAPoint( -m_dMidPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
      if ( m_nColorMethod == CM_LinearOpaque )
        m_lut->AddRGBAPoint( -m_dMinPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
      else
        m_lut->AddRGBAPoint( -m_dMinPoint, 0.4, 0.4, 0.4, 1 );
    }
    else
    {
      m_lut->AddRGBAPoint( -m_dMaxPoint, 1-m_colorMax[0]/255.0, m_colorMax[1]/255.0, 1-m_colorMax[2]/255.0, 1 );
      m_lut->AddRGBAPoint( -m_dMidPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
      if ( m_nColorMethod == CM_LinearOpaque )
        m_lut->AddRGBAPoint( -m_dMinPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
      else
        m_lut->AddRGBAPoint( -m_dMinPoint, 0.4, 0.4, 0.4, 1 );
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
        m_lut->AddRGBAPoint( m_dMinPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
      else
        m_lut->AddRGBAPoint( m_dMinPoint, 0.4, 0.4, 0.4, 1 );
      m_lut->AddRGBAPoint( m_dMidPoint, 1-m_colorMin[0]/255.0, m_colorMin[1]/255.0, 1-m_colorMin[2]/255.0, 1 );
      m_lut->AddRGBAPoint( m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
    }
    else
    {
      if ( m_nColorMethod == CM_LinearOpaque )
        m_lut->AddRGBAPoint( m_dMinPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
      else
        m_lut->AddRGBAPoint( m_dMinPoint, 0.4, 0.4, 0.4, 1 );
      m_lut->AddRGBAPoint( m_dMidPoint, m_colorMin[0]/255.0, m_colorMin[1]/255.0, m_colorMin[2]/255.0, 1 );
      m_lut->AddRGBAPoint( m_dMaxPoint, m_colorMax[0]/255.0, m_colorMax[1]/255.0, m_colorMax[2]/255.0, 1 );
    }
  }
 
  m_lut->Build();
}

int SurfaceOverlayProperties::GetColorMethod()
{
  return m_nColorMethod;
}

void SurfaceOverlayProperties::SetColorMethod( int nTh )
{
  m_nColorMethod = nTh;
}

void SurfaceOverlayProperties::SetMinPoint( double dValue )
{
  if ( dValue != m_dMinPoint )
  {
    m_dMinPoint = dValue;
  }
}
  
double SurfaceOverlayProperties::GetMinPoint()
{
  return m_dMinPoint;
}
  
  
void SurfaceOverlayProperties::SetMidPoint( double dValue )
{
  if ( dValue != m_dMidPoint )
  {
    m_dMidPoint = dValue;
  }
}
  
double SurfaceOverlayProperties::GetMidPoint()
{
  return m_dMidPoint;
}
  
void SurfaceOverlayProperties::SetMaxPoint( double dValue )
{
  if ( dValue != m_dMaxPoint )
  {
    m_dMaxPoint = dValue;
  }
}
  
double SurfaceOverlayProperties::GetMaxPoint()
{
  return m_dMaxPoint;
}

bool SurfaceOverlayProperties::GetColorInverse()
{
  return m_bColorInverse;
}

void SurfaceOverlayProperties::SetColorInverse( bool bInverse )
{
  m_bColorInverse = bInverse;
  SetColorScale( m_nColorScale );
}

bool SurfaceOverlayProperties::GetColorTruncate()
{
  return m_bColorTruncate;
}

void SurfaceOverlayProperties::SetColorTruncate( bool bTruncate )
{
  m_bColorTruncate = bTruncate;
  SetColorScale( m_nColorScale );
}


void SurfaceOverlayProperties::MapOverlayColor( unsigned char* colordata, int nPoints )
{
  MapOverlayColor( m_overlay->GetData(), colordata, nPoints );
}

void SurfaceOverlayProperties::MapOverlayColor( float* data, unsigned char* colordata, int nPoints )
{
  double c[3];
  if ( m_nColorMethod== CM_LinearOpaque )
  {
    for ( int i = 0; i < nPoints; i++ )
    {
      // map positive values
      if ( data[i] >= m_dMinPoint && !( m_bColorInverse && m_bColorTruncate ) )
      {
        double r = ( m_dMaxPoint - data[i] ) / ( m_dMaxPoint - m_dMidPoint );
        if ( r < 0 )
          r = 0;
        else if ( r > 1 )
          r = 1;
       
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
        double r = ( data[i] + m_dMaxPoint ) / ( m_dMaxPoint - m_dMidPoint ); 
        if ( r < 0 )
          r = 0;
        else if ( r > 1 )
          r = 1;
              
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
        if ( data[i] < m_dMidPoint )
        {
          double r = ( m_dMidPoint - data[i] ) / ( m_dMidPoint - m_dMinPoint );
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
          double r = ( m_dMaxPoint - data[i] ) / ( m_dMaxPoint - m_dMidPoint );
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
        if ( data[i] >= -m_dMidPoint )
        {
          double r = ( m_dMidPoint + data[i] ) / ( m_dMidPoint - m_dMinPoint );
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
          double r = ( m_dMaxPoint + data[i] ) / ( m_dMaxPoint - m_dMidPoint );
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

