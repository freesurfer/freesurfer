/**
 * @file  DialogGradientVolume.h
 * @brief Dialog to create gradient volume.
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



#include "DialogGradientVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "LayerMRI.h"
#include "VolumeFilterGradient.h"

BEGIN_EVENT_TABLE( DialogGradientVolume, wxDialog )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_APPLY" ),     DialogGradientVolume::OnApply )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CANCEL" ),    DialogGradientVolume::OnClose )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_SMOOTH" ),     DialogGradientVolume::OnCheckSmooth )
END_EVENT_TABLE()


DialogGradientVolume::DialogGradientVolume( LayerMRI* input, LayerMRI* output, wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_GRADIENT_VOLUME") );
  m_textSD      =   XRCCTRL( *this, "ID_TEXT_SD",       wxTextCtrl );
  m_checkSmooth =   XRCCTRL( *this, "ID_CHECK_SMOOTH",  wxCheckBox );
  m_btnApply    =   XRCCTRL( *this, "ID_BUTTON_APPLY",  wxButton );
  
  m_filter = new VolumeFilterGradient( input, output );
}

DialogGradientVolume::~DialogGradientVolume()
{
  delete m_filter;
}

void DialogGradientVolume::SetVolumes( LayerMRI* input, LayerMRI* output )
{
  m_filter->SetInputOutputVolumes( input, output );
}


void DialogGradientVolume::OnApply( wxCommandEvent& event )
{
  double dSD;
  if ( !m_textSD->GetValue().ToDouble( &dSD ) || dSD <= 0 )
  {
    wxMessageDialog dlg
      ( this, 
	_("Standard deviation must be valid and greater than 0."), 
	_("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  
  m_btnApply->Enable( false );
  m_filter->SetSmoothing( m_checkSmooth->IsChecked() );
  m_filter->SetStandardDeviation( dSD );
  m_filter->Update();
  
  m_btnApply->Enable( true );
  
  return;
}

void DialogGradientVolume::OnCheckSmooth( wxCommandEvent& event )
{
  m_textSD->Enable( event.IsChecked() );
}
