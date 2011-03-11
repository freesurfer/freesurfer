/**
 * @file  DialogNewWayPoints.h
 * @brief Dialog to create new way points.
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



#include "DialogNewWayPoints.h"
#include <wx/xrc/xmlres.h>
#include "stdlib.h"
#include "stdio.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "LayerPropertiesWayPoints.h"

BEGIN_EVENT_TABLE( DialogNewWayPoints, wxDialog )
  EVT_BUTTON     ( wxID_OK,         DialogNewWayPoints::OnOK )
END_EVENT_TABLE()


DialogNewWayPoints::DialogNewWayPoints( wxWindow* parent, 
					LayerCollection* col_mri )
{
  wxXmlResource::Get()->LoadDialog( this, parent, 
				    _("ID_DIALOG_NEW_WAYPOINTS") );
  m_textName        = XRCCTRL( *this, "ID_TEXT_NAME", wxTextCtrl );
  m_choiceTemplate  = XRCCTRL( *this, "ID_CHOICE_TEMPLATE", wxChoice );
  m_radioControlPoints  = XRCCTRL( *this, "ID_RADIO_CONTROL_POINTS", wxRadioButton );
  m_radioWayPoints      = XRCCTRL( *this, "ID_RADIO_WAY_POINTS", wxRadioButton );

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

DialogNewWayPoints::~DialogNewWayPoints()
{}

wxString DialogNewWayPoints::GetWayPointsName()
{
  return m_textName->GetValue().Trim( true ).Trim( false );
}

void DialogNewWayPoints::SetWayPointsName( const wxString& name )
{
  m_textName->SetValue( name );
}

LayerMRI* DialogNewWayPoints::GetTemplate()
{
  return ( LayerMRI* )( void* )m_choiceTemplate->
    GetClientData( m_choiceTemplate->GetSelection() );
}

void DialogNewWayPoints::OnOK( wxCommandEvent& event )
{
  if ( GetWayPointsName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Way points name cannot be empty."), 
			 _("Error"), 
			 wxOK );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}

int DialogNewWayPoints::GetType()
{
  if ( m_radioControlPoints->GetValue() )
    return LayerPropertiesWayPoints::ControlPoints;
  else
    return LayerPropertiesWayPoints::WayPoints;
}
