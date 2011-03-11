/**
 * @file  DialogLoadROI.h
 * @brief Dialog to load ROI/label data.
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



#include "DialogLoadROI.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "MyUtils.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LayerCollection.h"

BEGIN_EVENT_TABLE( DialogLoadROI, wxDialog )
  EVT_BUTTON    ( wxID_OK,                        DialogLoadROI::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_OPEN" ),      DialogLoadROI::OnButtonOpen )
END_EVENT_TABLE()


DialogLoadROI::DialogLoadROI( wxWindow* parent, bool bEnableResample )
{
  wxXmlResource::Get()->LoadDialog( this, parent, 
				    wxT("ID_DIALOG_LOAD_ROI") );
  m_choiceTemplate  = XRCCTRL( *this, "ID_CHOICE_TEMPLATE", wxChoice );
  m_textFileName    = XRCCTRL( *this, "ID_TEXT_FILENAME",   wxTextCtrl );
  
  LayerCollection* col_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
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
}

DialogLoadROI::~DialogLoadROI()
{}

wxArrayString DialogLoadROI::GetFileNames()
{
  return MyUtils::SplitString( m_textFileName->GetValue(), _(";") );
}


void DialogLoadROI::OnOK( wxCommandEvent& event )
{
  if ( GetFileNames().IsEmpty())
  {
    wxMessageDialog dlg( this, 
			 _("ROI file names cannot be empty."), 
			 _("Error"), 
			 wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}

wxString DialogLoadROI::GetTemplate()
{
  return m_choiceTemplate->GetStringSelection();
}

void DialogLoadROI::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg
    ( this, 
      _("Open ROI files"), 
      m_strLastDir, _(""),
      _("ROI files (*.label)|*.label|All files (*.*)|*.*"),
      wxFD_OPEN | wxFD_MULTIPLE );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxArrayString fns;
    dlg.GetPaths( fns );
    wxString text;
    for ( size_t i = 0; i < fns.GetCount(); i++ )
    {
      text += fns[i];
      if ( i != fns.GetCount()-1 )
        text += _(";");
    }
    m_textFileName->SetValue( text );
    m_textFileName->SetInsertionPointEnd();
 //   m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

