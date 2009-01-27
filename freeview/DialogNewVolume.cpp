/**
 * @file  DialogNewVolume.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:43:47 $
 *    $Revision: 1.2.2.2 $
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



#include "DialogNewVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/msgdlg.h>
#include "stdlib.h"
#include "stdio.h"
#include "LayerMRI.h"
#include "LayerCollection.h"

BEGIN_EVENT_TABLE( DialogNewVolume, wxDialog )
EVT_BUTTON     ( wxID_OK,         DialogNewVolume::OnOK )
EVT_TEXT_ENTER    ( XRCID( wxT( "ID_TEXT_NAME" ) ),  DialogNewVolume::OnTextEnter )
END_EVENT_TABLE()


DialogNewVolume::DialogNewVolume( wxWindow* parent, LayerCollection* col_mri )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_NEW_VOLUME") );
  m_checkCopyVoxel = XRCCTRL( *this, "ID_CHECK_COPY_VOXEL", wxCheckBox );
  m_textName = XRCCTRL( *this, "ID_TEXT_NAME", wxTextCtrl );
  m_choiceTemplate = XRCCTRL( *this, "ID_CHOICE_TEMPLATE", wxChoice );

  std::vector<Layer*> layers = col_mri->GetLayers();
  int nSel = 0;
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    m_choiceTemplate->Insert( layers[i]->GetName(), 0, (void*)layers[i] );
    if ( layers[i] == col_mri->GetActiveLayer() )
      nSel = i;
  }
  m_choiceTemplate->SetSelection( nSel );
  m_textName->SetFocus();
}

DialogNewVolume::~DialogNewVolume()
{}

wxString DialogNewVolume::GetVolumeName()
{
  return m_textName->GetValue().Trim( true ).Trim( false );
}

void DialogNewVolume::SetVolumeName( const wxString& name )
{
  m_textName->SetValue( name );
}

bool DialogNewVolume::GetCopyVoxel()
{
  return m_checkCopyVoxel->IsChecked();
}

void DialogNewVolume::SetCopyVoxel( bool bVoxel )
{
  m_checkCopyVoxel->SetValue( bVoxel );
}

LayerMRI* DialogNewVolume::GetTemplate()
{
  return ( LayerMRI* )( void* )m_choiceTemplate->GetClientData( m_choiceTemplate->GetSelection() );
}

void DialogNewVolume::OnOK( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, "Volume name can not be empty.", "Error", wxOK );
    dlg.ShowModal();
    return;
  }
  event.Skip();
}


void DialogNewVolume::OnTextEnter( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, "Volume name can not be empty.", "Error", wxOK );
    dlg.ShowModal();
    return;
  }
  EndModal( wxID_OK );
}

