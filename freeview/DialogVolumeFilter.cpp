/**
 * @file  DialogVolumeFilter.h
 * @brief Dialog to create gradient volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/26 19:04:05 $
 *    $Revision: 1.1 $
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



#include "DialogVolumeFilter.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include "LayerMRI.h"
#include "VolumeFilter.h"

BEGIN_EVENT_TABLE( DialogVolumeFilter, wxDialog )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_APPLY" ),     DialogVolumeFilter::OnApply )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CANCEL" ),    DialogVolumeFilter::OnClose )
END_EVENT_TABLE()


DialogVolumeFilter::DialogVolumeFilter( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_VOLUME_FILTER") );
  m_btnApply        =   XRCCTRL( *this, "ID_BUTTON_APPLY",  wxButton );
  m_spinKernelSize  =   XRCCTRL( *this, "ID_SPIN_KERNEL_SIZE",  wxSpinCtrl );
  m_textSigma       =   XRCCTRL( *this, "ID_TEXT_SIGMA",  wxTextCtrl );
}

DialogVolumeFilter::~DialogVolumeFilter()
{
}

void DialogVolumeFilter::SetFilter( VolumeFilter* filter )
{
  m_filter = filter;
  if ( filter )
  {
    m_spinKernelSize->SetValue( filter->GetKernelSize() );
    SetTitle( wxString("Apply ") + filter->GetName().c_str() + _(" Filter") );
  }
}

int DialogVolumeFilter::GetKernelSize()
{
  return m_spinKernelSize->GetValue();
}

void DialogVolumeFilter::SetSigma( double dvalue )
{
  m_textSigma->SetValue( ( wxString() << dvalue ) );
}

double DialogVolumeFilter::GetSigma()
{
  double dvalue = 0;
  m_textSigma->GetValue().ToDouble( &dvalue );
  return dvalue;
}

void DialogVolumeFilter::ShowSigma( bool bshow )
{
  m_textSigma->Show( bshow );
  XRCCTRL( *this, "ID_STATIC_SIGMA",  wxStaticText )->Show( bshow );
}

void DialogVolumeFilter::OnApply( wxCommandEvent& event )
{
  if ( GetKernelSize() <= 0 )
  {
    wxMessageDialog dlg( this, _("Kernel size must be greater than 0."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  else if ( m_textSigma->IsShown() && GetSigma() <= 0 )
  {
    wxMessageDialog dlg( this, _("Sigma must be greater than 0."), 
                         _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  EndModal( wxID_OK );
}

