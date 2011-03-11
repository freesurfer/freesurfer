/**
 * @file  WindowConnectivityConfiguration.h
 * @brief Overlay configuration window.
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
#include "WindowConnectivityConfiguration.h"
#include <wx/xrc/xmlres.h>
#include <wx/config.h>
#include "wxHistogramWidget.h"
#include "Listener.h"
#include "MainWindow.h"
#include "ConnectivityData.h"

#define RADIUS_MIN    0.1
#define RADIUS_MAX    10.

BEGIN_EVENT_TABLE( WindowConnectivityConfiguration, wxFrame )
	EVT_CLOSE           ( WindowConnectivityConfiguration::OnClose )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_DISPLAY" ),        WindowConnectivityConfiguration::OnChoiceDisplay )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_RADIUS" ),         WindowConnectivityConfiguration::OnChoiceRadius )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_RADIUS_MIN" ),     WindowConnectivityConfiguration::OnSliderRadiusMin )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_RADIUS_MAX" ),     WindowConnectivityConfiguration::OnSliderRadiusMax )
  EVT_TEXT            ( XRCID( "ID_TEXT_RADIUS_MIN" ),       WindowConnectivityConfiguration::OnTextRadiusMin )
  EVT_TEXT            ( XRCID( "ID_TEXT_RADIUS_MAX" ),       WindowConnectivityConfiguration::OnTextRadiusMax )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_THRESHOLD_MIN" ),  WindowConnectivityConfiguration::OnSliderThresholdMin )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_THRESHOLD_MAX" ),  WindowConnectivityConfiguration::OnSliderThresholdMax )
  EVT_TEXT            ( XRCID( "ID_TEXT_THRESHOLD_MIN" ),    WindowConnectivityConfiguration::OnTextThresholdMin )
  EVT_TEXT            ( XRCID( "ID_TEXT_THRESHOLD_MIN" ),    WindowConnectivityConfiguration::OnTextThresholdMin )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_ADD_ON" ),       WindowConnectivityConfiguration::OnCheckAddOn )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_EXPORT" ),         WindowConnectivityConfiguration::OnButtonExport )
END_EVENT_TABLE()

WindowConnectivityConfiguration::WindowConnectivityConfiguration( wxWindow* parent ) : Listener( "WindowConnectivityConfiguration" )
{
  m_conn = MainWindow::GetMainWindowPointer()->GetConnectivityData();
  
	wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_FRAME_CONNECTIVITY") );

  m_choiceDisplay       = XRCCTRL( *this, "ID_CHOICE_DISPLAY",            wxChoice );
  m_choiceRadius        = XRCCTRL( *this, "ID_CHOICE_RADIUS",             wxChoice );
  m_sliderRadiusMin     = XRCCTRL( *this, "ID_SLIDER_RADIUS_MIN",         wxSlider );
  m_sliderRadiusMax     = XRCCTRL( *this, "ID_SLIDER_RADIUS_MAX",         wxSlider );
  m_textRadiusMin       = XRCCTRL( *this, "ID_TEXT_RADIUS_MIN",           wxTextCtrl );
  m_textRadiusMax       = XRCCTRL( *this, "ID_TEXT_RADIUS_MAX",           wxTextCtrl );
  m_sliderThresholdMin  = XRCCTRL( *this, "ID_SLIDER_THRESHOLD_MIN",      wxSlider );
  m_sliderThresholdMax  = XRCCTRL( *this, "ID_SLIDER_THRESHOLD_MAX",      wxSlider );
  m_textThresholdMin    = XRCCTRL( *this, "ID_TEXT_THRESHOLD_MIN",        wxTextCtrl );
  m_textThresholdMax    = XRCCTRL( *this, "ID_TEXT_THRESHOLD_MAX",        wxTextCtrl );
  m_checkAddOn          = XRCCTRL( *this, "ID_CHECKBOX_ADD_ON",           wxCheckBox );
  m_btnExport           = XRCCTRL( *this, "ID_BUTTON_EXPORT",             wxButton );
 
  m_widgetsRadiusMax.push_back( m_sliderRadiusMax );
  m_widgetsRadiusMax.push_back( m_sliderRadiusMax );
  m_widgetsRadiusMax.push_back( XRCCTRL( *this, "ID_STATIC_RADIUS_MAX",        wxStaticText ) );

	wxConfigBase* config = wxConfigBase::Get();
	if ( config )
	{
		int x = config->Read( _T("/ConnectivityConfigurationWindow/PosX"), 280L );
    int y = config->Read( _T("/ConnectivityConfigurationWindow/PosY"), 30L );
		int x1 = 0, y1 = 0;
		if ( parent )
			parent->GetPosition( &x1, &y1 );
		Move( x1 + x, y1 + y );
	}
	
	UpdateUI();
}

void WindowConnectivityConfiguration::OnClose( wxCloseEvent& event )
{
	wxConfigBase* config = wxConfigBase::Get();
	if ( config && !IsIconized() )
	{
		int x, x2, y, y2, w, h;
		GetParent()->GetPosition( &x2, &y2 );
		GetPosition( &x, &y );
		GetSize( &w, &h );
    config->Write( _T("/ConnectivityConfigurationWindow/PosX"), (long) x - x2 );
    config->Write( _T("/ConnectivityConfigurationWindow/PosY"), (long) y - y2 );
	}

	Hide();
}


void WindowConnectivityConfiguration::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
	if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || 
        iMsg == "ActiveLayerChanged" || iMsg == "ActiveOverlayChanged" )
	{
    
  }
}

void WindowConnectivityConfiguration::ShowWindow( )
{
  UpdateUI();
	Show();
}


void WindowConnectivityConfiguration::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxFrame::OnInternalIdle();
}


void WindowConnectivityConfiguration::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void WindowConnectivityConfiguration::DoUpdateUI()
{
  if ( m_conn && m_conn->IsValid() )		
	{
    m_choiceDisplay->SetSelection( m_conn->GetDisplayMode() );
    m_choiceRadius->SetSelection( m_conn->GetVaryRadius()?1:0 );
    double range[2];
    m_conn->GetRadiusRange( range );
    m_sliderRadiusMin->SetValue( (int)( (range[0]-RADIUS_MIN)*100/(RADIUS_MAX-RADIUS_MIN) ) );
    UpdateTextValue( m_textRadiusMin, range[0] );
    m_sliderRadiusMax->SetValue( (int)( (range[1]-RADIUS_MIN)*100/(RADIUS_MAX-RADIUS_MIN) ) );
    UpdateTextValue( m_textRadiusMax, range[1] );
    
    m_sliderRadiusMax->Enable( m_conn->GetVaryRadius() );
    m_textRadiusMax->Enable( m_conn->GetVaryRadius() );
       
    double data_range[2];
    m_conn->GetDataRange( data_range );
    m_conn->GetThresholdRange( range );
    m_sliderThresholdMin->SetValue( (int)( (range[0]-data_range[0])*100/(data_range[1]-data_range[0]) ) );
    UpdateTextValue( m_textThresholdMin, range[0] );
    m_sliderThresholdMax->SetValue( (int)( (range[1]-data_range[0])*100/(data_range[1]-data_range[0]) ) );
    UpdateTextValue( m_textThresholdMax, range[1] );    
    
    m_checkAddOn->SetValue( m_conn->GetIncrementalDisplay() );
    m_checkAddOn->Enable( m_conn->GetDisplayMode() != ConnectivityData::DM_All );
	}
	else
	{
		Close();
	}
}

void WindowConnectivityConfiguration::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  wxString value_strg = ( wxString() << dvalue );
  if ( value_strg != ctrl->GetValue() && ( value_strg + _(".") ) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void WindowConnectivityConfiguration::OnChoiceDisplay( wxCommandEvent& event )
{
  if ( m_conn->IsValid() )
  {
    m_conn->SetDisplayMode( event.GetSelection() );
  }
}

void WindowConnectivityConfiguration::OnChoiceRadius( wxCommandEvent& event )
{
  if ( m_conn->IsValid() )
  {
    m_conn->SetVaryRadius( event.GetSelection() > 0 );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnSliderRadiusMin( wxScrollEvent& event )
{
  if ( m_conn->IsValid() )   
  {
    m_conn->SetRadius( m_sliderRadiusMin->GetValue()*(RADIUS_MAX-RADIUS_MIN)/100. + RADIUS_MIN );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnSliderRadiusMax( wxScrollEvent& event )
{
  if ( m_conn->IsValid() )   
  {
    m_conn->SetRadiusMax( m_sliderRadiusMax->GetValue()*(RADIUS_MAX-RADIUS_MIN)/100. + RADIUS_MIN );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnTextRadiusMin( wxCommandEvent& event )
{
  double dvalue;
  if ( m_conn->IsValid() && m_textRadiusMin->GetValue().ToDouble( &dvalue ) )   
  {
    m_conn->SetRadius( dvalue );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnTextRadiusMax( wxCommandEvent& event )
{
  double dvalue;
  if ( m_conn->IsValid() && m_textRadiusMax->GetValue().ToDouble( &dvalue ) )   
  {
    m_conn->SetRadiusMax( dvalue );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnSliderThresholdMin( wxScrollEvent& event )
{
  if ( m_conn->IsValid() )
  {
    double data_range[2];
    m_conn->GetDataRange( data_range );
    m_conn->SetThresholdMin( m_sliderThresholdMin->GetValue()*(data_range[1]-data_range[0])/100. + data_range[0] );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnSliderThresholdMax( wxScrollEvent& event )
{
  if ( m_conn->IsValid() )   
  {
    double data_range[2];
    m_conn->GetDataRange( data_range );
    m_conn->SetThresholdMax( m_sliderThresholdMax->GetValue()*(data_range[1]-data_range[0])/100. + data_range[0] );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnTextThresholdMin( wxCommandEvent& event )
{
  double dvalue;
  if ( m_conn->IsValid() && m_textThresholdMin->GetValue().ToDouble( &dvalue ) )   
  {
    m_conn->SetThresholdMin( dvalue );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnTextThresholdMax( wxCommandEvent& event )
{
  double dvalue;
  if ( m_conn->IsValid() && m_textThresholdMax->GetValue().ToDouble( &dvalue ) )   
  {
    m_conn->SetThresholdMax( dvalue );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnCheckAddOn( wxCommandEvent& event )
{
  if ( m_conn->IsValid() )
  {
    m_conn->SetIncrementalDisplay( m_checkAddOn->GetValue() );
    UpdateUI();
  }
}

void WindowConnectivityConfiguration::OnButtonExport( wxCommandEvent& event )
{
  if ( m_conn->IsValid() )
  {
    if ( !m_conn->HasAnySeeds() )
    {
      wxMessageDialog msg( this, _("There is no connectivity data to export."), 
                           _("Warning"), wxOK );
      msg.ShowModal();
      return;
    }
    
    wxFileDialog dlg( this, _("Export to file"), _(""), _(""),
                      _("GraphViz data files (*.txt)|*.txt"),
                      wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
    if ( dlg.ShowModal() == wxID_OK )
    {
      if ( !m_conn->Export( dlg.GetPath().c_str() ) )
      {
        wxMessageDialog msg( this, _("Can write to file."), 
                             _("Error"), wxOK );
        msg.ShowModal();
      }
    }    
  }
}
