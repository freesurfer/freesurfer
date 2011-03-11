/**
 * @file  wxHistogramWidget.cpp
 * @brief Widget to display histogram data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
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

#include <iostream>
#include <wx/wx.h>
#include "wxHistogramWidget.h"
#include <wx/dcmemory.h>
#include <wx/bitmap.h>
#include <math.h>
#include <stdio.h>
#include "MyUtils.h"


IMPLEMENT_DYNAMIC_CLASS( wxHistogramWidget, wxPanel )

BEGIN_EVENT_TABLE   ( wxHistogramWidget, wxPanel )
  EVT_PAINT         ( wxHistogramWidget::OnPaint )
  EVT_LEFT_DOWN     ( wxHistogramWidget::OnMouseButtonDown )
  EVT_RIGHT_DOWN    ( wxHistogramWidget::OnMouseButtonDown )
  EVT_MIDDLE_DOWN   ( wxHistogramWidget::OnMouseButtonDown )
END_EVENT_TABLE()

wxHistogramWidget::wxHistogramWidget() :
                   Broadcaster( "wxHistogramWidget" )
{
  Initialize();
}

wxHistogramWidget::wxHistogramWidget( wxWindow *parent, 
                                      wxWindowID id,
                                      const wxPoint& pos,
                                      const wxSize& size, 
                                      long style,
                                      const wxString& name ) :
                   Broadcaster( "wxHistogramWidget" )
{
  this->Create( parent, id, pos, size, style, name );
}


wxHistogramWidget::~wxHistogramWidget()
{
  if ( m_dInputData )
    delete[] m_dInputData;

  if ( m_nOutputData )
    delete[] m_nOutputData;
  
  if ( m_nColorTable )
    delete[] m_nColorTable;
}

bool wxHistogramWidget::Create( wxWindow *parent, 
								wxWindowID id,
								const wxPoint& pos,
								const wxSize& size, 
								long style,
								const wxString& name )
{
	if ( !wxPanel::Create( parent,id, pos, size, style | wxFULL_REPAINT_ON_RESIZE, name ) )
		return false;
	
  Initialize();
  
	return true;
}

void wxHistogramWidget::Initialize()
{
  m_dInputData = NULL;
  m_nInputSize = 0;
  m_nOutputData = NULL;
  m_nOutputSize = 0;
  m_bAutoRange = true;
  m_nNumberOfBins = 100;
  m_colorBackground = *wxWHITE;
  m_colorForeground = wxColour( 192, 192, 192 );
  m_nMaxCount = 0;
  m_nColorTable = NULL;
  
  m_rectGraph = wxRect( 50, 20, 100, 100 );
}

void wxHistogramWidget::GetOutputRange( double* dRange )
{
  dRange[0] = m_dOutputRange[0];
  dRange[1] = m_dOutputRange[1];
}

void wxHistogramWidget::SetOutputRange( double* dRange )
{
  m_dOutputRange[0] = dRange[0];
  m_dOutputRange[1] = dRange[1];
  UpdateData();
}

bool wxHistogramWidget::GetAutoRange()
{
  return m_bAutoRange;
}

void wxHistogramWidget::SetAutoRange( bool bRange )
{
  m_bAutoRange = bRange;
  if ( bRange )
  {
    m_dOutputRange[0] = m_dInputRange[0];
    m_dOutputRange[1] = m_dInputRange[1];
    UpdateData();
  }
}

int wxHistogramWidget::GetNumberOfBins()
{
  return m_nNumberOfBins;
}

void wxHistogramWidget::SetNumberOfBins( int nBins )
{
  m_nNumberOfBins = nBins;
  UpdateData();
}

int wxHistogramWidget::GetOutputSize()
{
  return m_nOutputSize;
}

void wxHistogramWidget::GetOutputData( int* buffer_out )
{
  memcpy( buffer_out, m_nOutputData, m_nOutputSize * sizeof( int ) );
}

void wxHistogramWidget::UpdateData( bool bRepaint )
{
	if ( m_dInputData == NULL || m_nNumberOfBins < 2 )
		return;
	
	m_dBinWidth = ( m_dOutputRange[1] - m_dOutputRange[0] ) / m_nNumberOfBins;
	if ( m_nOutputData )
		delete[] m_nOutputData;
	
	m_nOutputData = new int[m_nNumberOfBins];
	if ( !m_nOutputData )
	{
		std::cerr << "Can not allocate memory." << std::endl;
		return;
	}
	
  // calculate histogram data
	memset( m_nOutputData, 0, m_nNumberOfBins * sizeof( int ) );
	for ( long i = 0; i < m_nInputSize; i++ )
	{
		int n = (int)( ( m_dInputData[i] - m_dOutputRange[0] ) / m_dBinWidth );
		if ( n >= 0 && n < m_nNumberOfBins )
    {
			m_nOutputData[n] ++;
    }
	}
	
  // find max and second max
  m_nMaxCount = 0;
  int nSecondMax = 0;
  for ( int i = 0; i < m_nNumberOfBins; i++ )
  {
    if ( m_nMaxCount < m_nOutputData[i] )
      m_nMaxCount  = m_nOutputData[i];
    else if ( nSecondMax < m_nOutputData[i] )
      nSecondMax = m_nOutputData[i];
  }
    
  if ( m_nMaxCount > nSecondMax * 5 )
    m_nMaxCount = nSecondMax;
    
  // allocate color table
  if ( m_nColorTable )
    delete[] m_nColorTable;
  
  m_nColorTable = new unsigned char[ m_nNumberOfBins*4 ];
  UpdateColorTable();
  
	if ( bRepaint )
		Refresh();
}

void wxHistogramWidget::UpdateColorTable()
{ 
  for ( int i = 0; i < m_nNumberOfBins; i++ )
  {
    m_nColorTable[i*4] = m_colorForeground.Red();
    m_nColorTable[i*4+1] = m_colorForeground.Green();
    m_nColorTable[i*4+2] = m_colorForeground.Blue();
    m_nColorTable[i*4+3] = 255;
  }
}

void wxHistogramWidget::OnPaint( wxPaintEvent &event )
{
	wxPaintDC dc( this );
        
  wxSize sz = GetClientSize();
	if ( m_nOutputData )
	{	
    wxCoord x, y;  
    wxCoord nOrigin[2] = { m_rectGraph.x, m_rectGraph.y };
    wxCoord nCavWidth = sz.GetWidth() - m_rectGraph.x - 20;
    wxCoord nCavHeight = sz.GetHeight() - m_rectGraph.y - 20;
    m_rectGraph.SetWidth( nCavWidth );
    m_rectGraph.SetHeight( nCavHeight );
    
    // draw background
    dc.SetBrush( wxBrush( m_colorBackground ) );
    dc.SetPen( *wxTRANSPARENT_PEN );
    dc.DrawRectangle( m_rectGraph );
    
    // draw y metrics
    int nMetricInterval = 25;
    double dMetricStep = ((double)m_nMaxCount) / ( nCavHeight / nMetricInterval );
    dMetricStep = MyUtils::RoundToGrid( dMetricStep );
    double dMetricStart = 0;
    y = ( wxCoord )( m_rectGraph.GetBottom() );
    dc.SetPen( *wxBLACK_PEN );
    while ( y > m_rectGraph.GetTop() && dMetricStep > 0 )
    {
      if( y < m_rectGraph.GetBottom() )
      {
        dc.DrawLine( m_rectGraph.GetLeft(), y, m_rectGraph.GetLeft()-4, y );
        
        wxPen oldPen = dc.GetPen();
        dc.SetPen( wxPen( wxColour( 192, 192, 192 ), 1, wxDOT ) );
        dc.DrawLine( m_rectGraph.GetLeft(), y, m_rectGraph.GetRight(), y );
        dc.SetPen( oldPen );
      }
      wxString value_strg = ( wxString() << dMetricStart );
      wxCoord tx, ty;
      dc.GetTextExtent( value_strg, &tx, &ty );
      dc.DrawText( value_strg, m_rectGraph.GetLeft() - tx - 4, y - ty/2 );
      
      dMetricStart += dMetricStep;
      y = ( wxCoord )( m_rectGraph.GetBottom() - dMetricStart / m_nMaxCount *nCavHeight );
    }
    
    // draw bars 
    double dStepWidth = ( (double) nCavWidth) / m_nNumberOfBins;
    x = nOrigin[0];
    for ( int i = 0; i < m_nNumberOfBins; i++ )
    {
      dc.SetPen( wxPen( wxColour( m_nColorTable[i*4],  m_nColorTable[i*4+1], m_nColorTable[i*4+2] ) ) );
      dc.SetBrush( wxBrush( wxColour( m_nColorTable[i*4],  m_nColorTable[i*4+1], m_nColorTable[i*4+2] ) ) );
      y = (wxCoord)( nOrigin[1] + nCavHeight * ( 1.0 - (double)m_nOutputData[i] / m_nMaxCount ) );
      wxCoord h = (wxCoord)( (double)m_nOutputData[i] / m_nMaxCount * nCavHeight );
      if ( y < nOrigin[1] )
      {
        y = nOrigin[1]; 
        h = nCavHeight;
      }
      wxCoord x2 = (wxCoord)( nOrigin[0] + dStepWidth * (i+1) );
      
      dc.DrawRectangle( x, y, x2-x, h );            
      x = x2;
    }
        
    // draw markers
    for ( size_t i = 0; i < m_markers.size(); i++ )
    {
      switch ( m_markers[i].style )
      {
        case LineMarker::Solid:
          dc.SetPen( wxPen( m_markers[i].color, 1, wxSOLID ) );
          break;
        case LineMarker::Dash:
          dc.SetPen( wxPen( m_markers[i].color, 1, wxSHORT_DASH ) );
          break;
        case LineMarker::Dot:
          dc.SetPen( wxPen( m_markers[i].color, 1, wxDOT ) );
          break;
      }
          
      x = ( wxCoord )( nOrigin[0] + nCavWidth * ( m_markers[i].position - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] ) );
      y = nOrigin[1];
      if ( x > nOrigin[0] && x < nOrigin[0] + nCavWidth ) 
        dc.DrawLine( x, y, x, y + nCavHeight );
    }
    
    // draw axis 
    dc.SetBrush( *wxTRANSPARENT_BRUSH );
    dc.SetPen( *wxBLACK_PEN );
    dc.DrawRectangle( m_rectGraph );
    
    // draw x metrics
    nMetricInterval = 50;
    dMetricStep = ( m_dOutputRange[1] - m_dOutputRange[0] ) / ( nCavWidth / nMetricInterval );
    dMetricStep = MyUtils::RoundToGrid( dMetricStep );
    dMetricStart = (int)( ( m_dOutputRange[0] / dMetricStep ) ) * dMetricStep;
    if ( m_dOutputRange[0] < 0 )
      dMetricStart -= dMetricStep;
    x = ( wxCoord )( nOrigin[0] + ( dMetricStart - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] ) *nCavWidth );
    while ( x < m_rectGraph.GetRight() && dMetricStep > 0 )
    {
      if( x >= m_rectGraph.GetLeft() )
      {
        dc.DrawLine( x, m_rectGraph.GetBottom(), x, m_rectGraph.GetBottom()+4 );
      } 
      wxString value_strg = ( wxString() << dMetricStart );
      wxCoord tx, ty;
      dc.GetTextExtent( value_strg, &tx, &ty );
      if ( x - tx / 2 > 0 )
        dc.DrawText( value_strg, x - tx / 2, m_rectGraph.GetBottom() + 5 );
      
      dMetricStart += dMetricStep;
      if ( fabs( dMetricStart ) < 1e-10 )
        dMetricStart = 0;
      
      x = ( wxCoord )(m_rectGraph.GetLeft() + ( dMetricStart - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] )*nCavWidth );
    }
	}
	else
	{
		event.Skip();
	}
}

void wxHistogramWidget::SetForegroundColor( const wxColour& color )
{
  m_colorForeground = color;
}

void wxHistogramWidget::SetColorTableData( unsigned char* colortable, bool bRefresh )
{
  memcpy( m_nColorTable, colortable, m_nNumberOfBins * 4 );
  if ( bRefresh )
    Refresh();
}

void wxHistogramWidget::SetMarkers( const MarkerArray& markers, bool bRefresh  )
{
  m_markers = markers;
  if ( bRefresh )
    Refresh();
}

void wxHistogramWidget::OnMouseButtonDown( wxMouseEvent& event )
{
  if ( m_rectGraph.Contains( event.GetX(), event.GetY() ) )
  {
    double dValue = ( event.GetX() - m_rectGraph.GetX() ) * ( m_dOutputRange[1] - m_dOutputRange[0] ) / m_rectGraph.GetWidth() + m_dOutputRange[0];
    if ( event.LeftDown() )
    {
      this->SendBroadcast( "LeftButtonPressed", &dValue, this );
    }
    else if ( event.MiddleDown() )
    {
      this->SendBroadcast( "MiddleButtonPressed", &dValue, this );
    }
    else if ( event.RightDown() )
    {
      this->SendBroadcast( "RightButtonPressed", &dValue, this );
    }
  }
    
  event.Skip();
}
