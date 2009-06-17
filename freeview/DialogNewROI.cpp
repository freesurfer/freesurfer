/**
 * @file  DialogNewROI.h
 * @brief Dialog to create new ROI.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/06/17 20:41:17 $
 *    $Revision: 1.10 $
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



#include "DialogNewROI.h"
#include <wx/xrc/xmlres.h>
#include "stdlib.h"
#include "stdio.h"
#include "LayerMRI.h"
#include "LayerCollection.h"

BEGIN_EVENT_TABLE( DialogNewROI, wxDialog )
  EVT_BUTTON     ( wxID_OK,         DialogNewROI::OnOK )
  EVT_TEXT_ENTER ( XRCID( "ID_TEXT_NAME" ),  DialogNewROI::OnTextEnter )
END_EVENT_TABLE()


DialogNewROI::DialogNewROI( wxWindow* parent, LayerCollection* col_mri )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_NEW_ROI") );
  m_textName        = XRCCTRL( *this, "ID_TEXT_NAME", wxTextCtrl );
  m_choiceTemplate  = XRCCTRL( *this, "ID_CHOICE_TEMPLATE", wxChoice );

  std::vector<Layer*> layers = col_mri->GetLayers();
  int nSel = 0;
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    m_choiceTemplate->Insert( wxString::FromAscii(layers[i]->GetName()), 
			      0, (void*)layers[i] );
    if ( layers[i] == col_mri->GetActiveLayer() )
      nSel = i;
  }
  m_choiceTemplate->SetSelection( nSel );
  m_textName->SetFocus();
}

DialogNewROI::~DialogNewROI()
{}

wxString DialogNewROI::GetROIName()
{
  return m_textName->GetValue().Trim( true ).Trim( false );
}

void DialogNewROI::SetROIName( const wxString& name )
{
  m_textName->SetValue( name );
}

LayerMRI* DialogNewROI::GetTemplate()
{
  return ( LayerMRI* )( void* )m_choiceTemplate->
    GetClientData( m_choiceTemplate->GetSelection() );
}

void DialogNewROI::OnOK( wxCommandEvent& event )
{
  if ( GetROIName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("ROI name cannot be empty."), 
			 _("Error"), 
			 wxOK );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}


void DialogNewROI::OnTextEnter( wxCommandEvent& event )
{
  if ( GetROIName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("ROI name cannot be empty."), 
			 _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  EndModal( wxID_OK );
}
