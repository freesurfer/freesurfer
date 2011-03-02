/**
 * @file  DialogNewVolume.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.12 $
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



#include "DialogNewVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/msgdlg.h>
#include "stdlib.h"
#include "stdio.h"
#include "LayerMRI.h"
#include "LayerCollection.h"

BEGIN_EVENT_TABLE( DialogNewVolume, wxDialog )
  EVT_BUTTON     ( wxID_OK,                 DialogNewVolume::OnOK )
  EVT_TEXT_ENTER ( XRCID( "ID_TEXT_NAME" ), DialogNewVolume::OnTextEnter )
END_EVENT_TABLE()


DialogNewVolume::DialogNewVolume( wxWindow* parent, LayerCollection* col_mri )
{
  wxXmlResource::Get()->LoadDialog
    ( this, parent, wxT("ID_DIALOG_NEW_VOLUME") );
  m_checkCopyVoxel  = XRCCTRL( *this, "ID_CHECK_COPY_VOXEL", wxCheckBox );
  m_textName        = XRCCTRL( *this, "ID_TEXT_NAME", wxTextCtrl );
  m_choiceTemplate  = XRCCTRL( *this, "ID_CHOICE_TEMPLATE", wxChoice );
  m_choiceDataType  = XRCCTRL( *this, "ID_CHOICE_DATA_TYPE", wxChoice );

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
  return ( LayerMRI* )( void* )m_choiceTemplate->
    GetClientData( m_choiceTemplate->GetSelection() );
}

int DialogNewVolume::GetDataType()
{
  if ( m_choiceDataType->GetSelection() == (int)m_choiceDataType->GetCount()-1 )
    return GetTemplate()->GetDataType();
  else
    return m_choiceDataType->GetSelection();
}

void DialogNewVolume::OnOK( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Volume name can not be empty."), 
			 _("Error"), 
			 wxOK );
    dlg.ShowModal();
    return;
  }
  event.Skip();
}


void DialogNewVolume::OnTextEnter( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Volume name can not be empty."), 
			 _("Error"), 
			 wxOK );
    dlg.ShowModal();
    return;
  }
  EndModal( wxID_OK );
}

