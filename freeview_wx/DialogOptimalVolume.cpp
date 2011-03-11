/**
 * @file  DialogOptimalVolume.h
 * @brief Dialog to compute optimal volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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



#include "DialogOptimalVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/msgdlg.h>
#include "stdlib.h"
#include "stdio.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LayerCollection.h"

BEGIN_EVENT_TABLE( DialogOptimalVolume, wxDialog )
  EVT_BUTTON      ( wxID_OK,                  DialogOptimalVolume::OnOK )
  EVT_TEXT_ENTER  ( XRCID( "ID_TEXT_NAME" ),  DialogOptimalVolume::OnTextEnter )
END_EVENT_TABLE()


DialogOptimalVolume::DialogOptimalVolume( wxWindow* parent, 
					  LayerCollection* col_mri )
{
  wxXmlResource::Get()->LoadDialog( this, parent, 
				    wxT("ID_DIALOG_OPTIMAL_VOLUME") );
  m_textName          = XRCCTRL( *this, "ID_TEXT_NAME", wxTextCtrl );
  m_choiceLabelVolume = XRCCTRL( *this, "ID_CHOICE_LABEL_VOLUME", wxChoice );
  m_listBoxLayers     = XRCCTRL( *this, "ID_LISTBOX_LAYERS", wxCheckListBox );

  std::vector<Layer*> layers = col_mri->GetLayers();
  int nSel = 0;
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    m_choiceLabelVolume->Append( wxString::FromAscii(layers[i]->GetName()), 
				 (void*)layers[i] );
    if ( layers[i] == col_mri->GetActiveLayer() )
      nSel = i;
  }
  m_choiceLabelVolume->SetSelection( nSel );
  m_textName->SetFocus();

  for ( int i = 0; i < col_mri->GetNumberOfLayers(); i++ )
  {
    LayerMRI* layer = (LayerMRI*)col_mri->GetLayer( i );
    m_listBoxLayers->Append( wxString::FromAscii(layer->GetName()), 
			     (void*)layer );
    m_listBoxLayers->Check
      ( i, 
	layer->GetProperties()->GetColorMap() != LayerPropertiesMRI::LUT );
  }
}

DialogOptimalVolume::~DialogOptimalVolume()
{}

wxString DialogOptimalVolume::GetVolumeName()
{
  return m_textName->GetValue().Trim( true ).Trim( false );
}

void DialogOptimalVolume::SetVolumeName( const wxString& name )
{
  m_textName->SetValue( name );
}

LayerMRI* DialogOptimalVolume::GetLabelVolume()
{
  return ( LayerMRI* )( void* )m_choiceLabelVolume->
    GetClientData( m_choiceLabelVolume->GetSelection() );
}

void DialogOptimalVolume::OnOK( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Volume name cannot be empty."), 
			 _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  event.Skip();
}


void DialogOptimalVolume::OnTextEnter( wxCommandEvent& event )
{
  if ( GetVolumeName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Volume name cannot be empty."), 
			 _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  else if ( GetSelectedLayers().size() < 2 )
  {
    wxMessageDialog dlg
      ( this, 
	_("Must select at least two volumes for calculation."), 
	_("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
// EndModal( wxID_OK );
}

std::vector<LayerMRI*> DialogOptimalVolume::GetSelectedLayers()
{
  std::vector<LayerMRI*> layers;
  for ( size_t i = 0; i < m_listBoxLayers->GetCount(); i++ )
  {
    if ( m_listBoxLayers->IsChecked( i ) )
    {
      LayerMRI* layer = 
	( LayerMRI* )( void* )m_listBoxLayers->GetClientData( i );
      if ( layer )
        layers.push_back( layer );
    }
  }
  return layers;
}


