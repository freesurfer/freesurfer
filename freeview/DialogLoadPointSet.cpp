/**
 * @file  DialogLoadPointSet.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/02/16 20:49:01 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008-2009,
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



#include "DialogLoadPointSet.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "LayerPropertiesWayPoints.h"

BEGIN_EVENT_TABLE( DialogLoadPointSet, wxDialog )
  EVT_BUTTON    ( wxID_OK, DialogLoadPointSet::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_OPEN" ),        DialogLoadPointSet::OnButtonOpen )
END_EVENT_TABLE()


DialogLoadPointSet::DialogLoadPointSet( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_LOAD_POINTSET") );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_radioWayPoints      = XRCCTRL( *this, "ID_RADIO_WAY_POINTS",      wxRadioButton );
  m_radioControlPoints  = XRCCTRL( *this, "ID_RADIO_CONTROL_POINTS",  wxRadioButton );
}

DialogLoadPointSet::~DialogLoadPointSet()
{}

wxString DialogLoadPointSet::GetFileName()
{
  return m_textFileName->GetValue().Trim( true ).Trim( false );
}


void DialogLoadPointSet::OnOK( wxCommandEvent& event )
{
  if ( GetFileName().IsEmpty() )
  {
    wxMessageDialog dlg
      ( this, 
	_("Point set file name can not be empty."), 
	_("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}

void DialogLoadPointSet::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg
    ( this, 
      _("Select point set file"), 
      m_strLastDir, 
      _(""),
      _("All files (*.*)|*.*"),
      wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textFileName->ChangeValue( dlg.GetPath() );
    m_textFileName->SetInsertionPointEnd();
    m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
    m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

int DialogLoadPointSet::GetPointSetType()
{
  if (m_radioControlPoints->GetValue() )
    return LayerPropertiesWayPoints::ControlPoints;
  else if ( m_radioWayPoints->GetValue() )
    return LayerPropertiesWayPoints::WayPoints;
  else
    return -1;
}
