/**
 * @file  wxColorIndicator.cpp
 * @brief wxColorIndicator class.
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


#include <wx/wx.h>
#include "wxColorIndicator.h"
#include <string>

IMPLEMENT_DYNAMIC_CLASS( wxColorIndicator, wxPanel )


BEGIN_EVENT_TABLE( wxColorIndicator, wxPanel )
EVT_PAINT  ( wxColorIndicator::OnPaint )
END_EVENT_TABLE()

wxColorIndicator::wxColorIndicator()
{}

wxColorIndicator::wxColorIndicator( wxWindow *parent,
                                    wxWindowID id,
                                    const wxPoint& pos,
                                    const wxSize& size,
                                    long style,
                                    const wxString& name )
{
  this->Create( parent, id, pos, size, style, name );
}


wxColorIndicator::~wxColorIndicator()
{}



bool wxColorIndicator:: Create( wxWindow *parent,
                                wxWindowID id,
                                const wxPoint& pos,
                                const wxSize& size,
                                long style,
                                const wxString& name )
{
  if ( !wxPanel::Create( parent,id, pos, size, style, name ) )
    return false;

  m_colorFill = wxNullColour;
  m_colorBorder = *wxBLACK;

  return true;
}

void wxColorIndicator::OnPaint( wxPaintEvent &event )
{
  wxPaintDC dc( this );

  if ( m_colorFill != wxNullColour )
  {
    wxSize sz = GetClientSize();
    dc.SetBrush( wxBrush( m_colorFill ) );
    dc.SetPen( wxPen( m_colorBorder ) );
    dc.DrawRoundedRectangle( 0, 0, sz.GetWidth(), sz.GetHeight(), 2 );
  }
  else
  {
    event.Skip();
  }
}

wxColour wxColorIndicator::GetColor()
{
  return m_colorFill;
}

void wxColorIndicator::SetColor( const wxColour& color )
{
  m_colorFill = color;
  Refresh();
}

wxColour wxColorIndicator::GetBorderColor()
{
  return m_colorBorder;
}

void wxColorIndicator::SetBorderColor( const wxColour& color )
{
  m_colorBorder = color;
  Refresh();
}

