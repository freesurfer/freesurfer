/**
 * @file  WindowOverlayConfiguration.h
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
#include "WindowOverlayConfiguration.h"
#include <wx/xrc/xmlres.h>
#include <wx/html/htmlwin.h>
#include <wx/config.h>
#include <wx/fs_mem.h>
#include "wxHistogramWidget.h"
#include "Listener.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerSurface.h"
#include "LayerPropertiesSurface.h"
#include "SurfaceOverlay.h"
#include "SurfaceOverlayProperties.h"
#include <vtkImageData.h>

BEGIN_EVENT_TABLE( WindowOverlayConfiguration, wxFrame )
  EVT_CLOSE          ( WindowOverlayConfiguration::OnClose )
  EVT_BUTTON         ( XRCID( "ID_BUTTON_CLOSE" ),          WindowOverlayConfiguration::OnButtonClose )
  EVT_BUTTON         ( XRCID( "ID_BUTTON_APPLY" ),          WindowOverlayConfiguration::OnButtonApply )
  EVT_CHOICE         ( XRCID( "ID_CHOICE_OVERLAY" ),        WindowOverlayConfiguration::OnChoiceOverlay )
  EVT_COMMAND_SCROLL ( XRCID( "ID_SLIDER_OPACITY" ),        WindowOverlayConfiguration::OnSliderOpacity )
  EVT_TEXT           ( XRCID( "ID_TEXT_OPACITY_CTRL" ),     WindowOverlayConfiguration::OnTextOpacity )
  EVT_TEXT           ( XRCID( "ID_TEXT_THRESHOLD_MIN" ),    WindowOverlayConfiguration::OnTextThresholdChanged )
  EVT_TEXT           ( XRCID( "ID_TEXT_THRESHOLD_MID" ),    WindowOverlayConfiguration::OnUpdateGraph )
  EVT_TEXT           ( XRCID( "ID_TEXT_THRESHOLD_MAX" ),    WindowOverlayConfiguration::OnTextThresholdChanged )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_GREEN_RED" ), WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_HEAT" ),      WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_BLUE_RED" ),  WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_COLOR_WHEEL" ),   WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_RYGB_WHEEL" ),    WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_COLOR_TWO_COND_GR" ),   WindowOverlayConfiguration::OnUpdateGraph )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_THRESHOLD_LINEAR" ),    WindowOverlayConfiguration::OnUpdateThreshold )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_THRESHOLD_LINEAR_OPAQUE" ),
                       WindowOverlayConfiguration::OnUpdateThreshold )
  EVT_RADIOBUTTON    ( XRCID( "ID_RADIO_THRESHOLD_PIECEWISE" ), WindowOverlayConfiguration::OnUpdateThreshold )
  EVT_CHECKBOX       ( XRCID( "ID_CHECKBOX_COLOR_INVERSE" ),    WindowOverlayConfiguration::OnUpdateGraph )
  EVT_CHECKBOX       ( XRCID( "ID_CHECKBOX_COLOR_TRUNCATE" ),   WindowOverlayConfiguration::OnUpdateGraph )
END_EVENT_TABLE()

WindowOverlayConfiguration::WindowOverlayConfiguration( wxWindow* parent ) : 
  Listener( "WindowOverlayConfiguration" ),
  m_textOpacity( NULL )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_FRAME_OVERLAY") );
//  m_histogramWidget  = XRCCTRL( *this, "ID_HISTOGRAM",                wxHistogramWidget );
  m_choiceOverlay    = XRCCTRL( *this, "ID_CHOICE_OVERLAY",           wxChoice );
  m_sliderOpacity    = XRCCTRL( *this, "ID_SLIDER_OPACITY",           wxSlider );
  m_textOpacity      = XRCCTRL( *this, "ID_TEXT_OPACITY_CTRL",        wxTextCtrl );
  
  m_radioGreenRed    = XRCCTRL( *this, "ID_RADIO_COLOR_GREEN_RED",    wxRadioButton );
  m_radioHeat        = XRCCTRL( *this, "ID_RADIO_COLOR_HEAT",         wxRadioButton );
  m_radioBlueRed     = XRCCTRL( *this, "ID_RADIO_COLOR_BLUE_RED",     wxRadioButton );
  m_radioColorWheel  = XRCCTRL( *this, "ID_RADIO_COLOR_COLOR_WHEEL",  wxRadioButton );
  m_radioRYGBWheel   = XRCCTRL( *this, "ID_RADIO_COLOR_RYGB_WHEEL",   wxRadioButton );
  m_radioTwoCondGR   = XRCCTRL( *this, "ID_RADIO_COLOR_TWO_COND_GR",  wxRadioButton );
  
  m_textThresholdMin = XRCCTRL( *this, "ID_TEXT_THRESHOLD_MIN",       wxTextCtrl );
  m_textThresholdMid = XRCCTRL( *this, "ID_TEXT_THRESHOLD_MID",       wxTextCtrl );
  m_textThresholdMax = XRCCTRL( *this, "ID_TEXT_THRESHOLD_MAX",       wxTextCtrl );
  
  m_radioThresholdLinear        = XRCCTRL( *this, "ID_RADIO_THRESHOLD_LINEAR",        wxRadioButton );
  m_radioThresholdLinearOpaque  = XRCCTRL( *this, "ID_RADIO_THRESHOLD_LINEAR_OPAQUE", wxRadioButton );
  m_radioThresholdPiecewise     = XRCCTRL( *this, "ID_RADIO_THRESHOLD_PIECEWISE",     wxRadioButton );
  
  m_checkColorInverse           = XRCCTRL( *this, "ID_CHECKBOX_COLOR_INVERSE",        wxCheckBox );
  m_checkColorTruncate          = XRCCTRL( *this, "ID_CHECKBOX_COLOR_TRUNCATE",       wxCheckBox );
 
  // workaround for wxWidgets 2.9.1 & wxformbuilder bug
  wxPanel* cpanel = XRCCTRL( *this, "ID_PANEL_HISTOGRAM", wxPanel );
  if ( cpanel )
  { 
    wxBoxSizer* sizer = new wxBoxSizer( wxHORIZONTAL );
    cpanel->SetSizer( sizer );
    m_histogramWidget = new wxHistogramWidget( cpanel );
    sizer->Add( m_histogramWidget, 1, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );
    m_histogramWidget->SetMinSize( wxSize( -1, 100 ) );
    m_histogramWidget->AddListener( this );
  }
      
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int x = config->Read( _T("/OverlayConfigurationWindow/PosX"), 280L );
    int y = config->Read( _T("/OverlayConfigurationWindow/PosY"), 30L );
    int w = config->Read( _T("/OverlayConfigurationWindow/Width"), 380L );
    int h = config->Read( _T("/OverlayConfigurationWindow/Height"), 520L );
    SetSize( w, h );
    int x1 = 0, y1 = 0;
    if ( parent )
      parent->GetPosition( &x1, &y1 );
    Move( x1 + x, y1 + y );
        
    m_histogramWidget->SetNumberOfBins( 200 ); //config->Read( _T("/OverlayConfigurationWindow/NumberOfBins"), 200L ) );
  }
  
  m_layerSurface = NULL;

  MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->AddListener( this );
  
  UpdateUI();
}

void WindowOverlayConfiguration::OnClose( wxCloseEvent& event )
{
  wxConfigBase* config = wxConfigBase::Get();
  if ( config && !IsIconized() )
  {
    int x, x2, y, y2, w, h;
    GetParent()->GetPosition( &x2, &y2 );
    GetPosition( &x, &y );
    GetSize( &w, &h );
    config->Write( _T("/OverlayConfigurationWindow/PosX"), (long) x - x2 );
    config->Write( _T("/OverlayConfigurationWindow/PosY"), (long) y - y2 );
    config->Write( _T("/OverlayConfigurationWindow/Width"), (long) w );
    config->Write( _T("/OverlayConfigurationWindow/Height"), (long) h );
  }

  Hide();
}


void WindowOverlayConfiguration::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || 
       iMsg == "ActiveLayerChanged" || iMsg == "ActiveOverlayChanged" )
  {
    m_layerSurface = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
    if ( m_layerSurface )
    {
      double* rgb = m_layerSurface->GetProperties()->GetBinaryColor();
      m_histogramWidget->SetForegroundColor( wxColour( (int)(rgb[0]*255), (int)(rgb[0]*255), (int)(rgb[0]*255) ) );
    }
    UpdateUI( true );   // force update
    UpdateGraph();
  }
  else if ( sender == m_histogramWidget )
  {
    double pos = fabs( *((double*)iData ) );
    if ( iMsg == "LeftButtonPressed" )
    {
      UpdateTextValue( m_textThresholdMin, pos );
    }
    else if ( iMsg == "MiddleButtonPressed" )
    {
      UpdateTextValue( m_textThresholdMid, pos );
    }
    else if ( iMsg == "RightButtonPressed" )
    {
      UpdateTextValue( m_textThresholdMax, pos );
    }
    UpdateThresholdChanges();
  }
}

void WindowOverlayConfiguration::ShowWindow( LayerSurface* default_layer )
{
  if ( m_layerSurface != default_layer )
  {
    m_layerSurface = default_layer;
  }
  
  UpdateUI();
  Show();
}


void WindowOverlayConfiguration::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxFrame::OnInternalIdle();
}


void WindowOverlayConfiguration::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void WindowOverlayConfiguration::DoUpdateUI()
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )   
  {
  //  m_choiceLayer->SetSelection( n );
    SurfaceOverlayProperties* p = m_layerSurface->GetActiveOverlay()->GetProperties();
    m_sliderOpacity   ->SetValue( (int)( p->GetOpacity() * 100 ) );
    UpdateTextValue( m_textOpacity, p->GetOpacity() );
    
    m_radioGreenRed   ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_GreenRed );
    m_radioHeat       ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_Heat );
    m_radioBlueRed    ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_BlueRed );
    m_radioColorWheel ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_ColorWheel );
    m_radioRYGBWheel  ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_RYGBWheel );
    m_radioTwoCondGR  ->SetValue( p->GetColorScale() == SurfaceOverlayProperties::CS_TwoCondGR );
    
    UpdateTextValue( m_textThresholdMin, p->GetMinPoint() );
    UpdateTextValue( m_textThresholdMax, p->GetMaxPoint() );
    UpdateTextValue( m_textThresholdMid, p->GetMidPoint() );
    
    m_radioThresholdLinear        ->SetValue( p->GetColorMethod() == SurfaceOverlayProperties::CM_Linear );
    m_radioThresholdLinearOpaque  ->SetValue( p->GetColorMethod() == SurfaceOverlayProperties::CM_LinearOpaque );
    m_radioThresholdPiecewise     ->SetValue( p->GetColorMethod() == SurfaceOverlayProperties::CM_Piecewise );
    
    m_checkColorInverse   ->SetValue( p->GetColorInverse() );
    m_checkColorTruncate  ->SetValue( p->GetColorTruncate() );
    
    m_textThresholdMid->Enable( m_radioThresholdPiecewise->GetValue() );
    XRCCTRL( *this, "ID_STATIC_MID_POINT", wxStaticText )->Enable( m_radioThresholdPiecewise->GetValue() );
  }
  else
  {
    Close();
  }
}

void WindowOverlayConfiguration::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  wxString value_strg = ( wxString() << dvalue );
  if ( value_strg != ctrl->GetValue() && ( value_strg + _(".") ) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void WindowOverlayConfiguration::OnChoiceOverlay( wxCommandEvent& event )
{
  if ( !m_layerSurface || !m_layerSurface->GetActiveOverlay() )
  {
    Close();
    return;
  }
//  m_activeOverlay = (Layer*)(void*)m_choiceLayer->GetClientData( event.GetSelection() );
  UpdateUI();
}

void WindowOverlayConfiguration::OnButtonApply( wxCommandEvent& event )
{
  if ( !m_layerSurface || !m_layerSurface->GetActiveOverlay() )
  {
    return;
  }
  
  SurfaceOverlayProperties* p = m_layerSurface->GetActiveOverlay()->GetProperties();
  if ( UpdateOverlayProperties( p ) )
    m_layerSurface->UpdateOverlay( true );  // ask to redraw
}

bool WindowOverlayConfiguration::UpdateOverlayProperties( SurfaceOverlayProperties* p )
{
  p->SetOpacity( m_sliderOpacity->GetValue() / 100.0 );
  
  wxString strg = m_textThresholdMin->GetValue();
  double dValue;
  if ( strg.ToDouble( &dValue ) )
    p->SetMinPoint( dValue );
  else
    return false;
  
  strg = m_textThresholdMid->GetValue();
  if ( strg.ToDouble( &dValue ) )
    p->SetMidPoint( dValue );
  else
    return false;
  
  strg = m_textThresholdMax->GetValue();
  if ( strg.ToDouble( &dValue ) )
    p->SetMaxPoint( dValue );
  else 
    return false;
  
  if ( m_radioGreenRed->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_GreenRed );
  else if ( m_radioHeat->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_Heat );
  else if ( m_radioBlueRed->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_BlueRed );
  else if ( m_radioColorWheel->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_ColorWheel );
  else if ( m_radioRYGBWheel->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_RYGBWheel );
  else if ( m_radioTwoCondGR->GetValue() )
    p->SetColorScale( SurfaceOverlayProperties::CS_TwoCondGR );
  
  if ( m_radioThresholdLinear->GetValue() )
    p->SetColorMethod( SurfaceOverlayProperties::CM_Linear );
  else if ( m_radioThresholdLinearOpaque->GetValue() )
    p->SetColorMethod( SurfaceOverlayProperties::CM_LinearOpaque );
  else if ( m_radioThresholdPiecewise->GetValue() )
    p->SetColorMethod( SurfaceOverlayProperties::CM_Piecewise );
  
  p->SetColorInverse( m_checkColorInverse->GetValue() );
  p->SetColorTruncate( m_checkColorTruncate->GetValue() );
  
  return true;
}

void WindowOverlayConfiguration::UpdateGraph()
{
  if ( m_layerSurface && m_layerSurface->GetActiveOverlay() )
  {
    SurfaceOverlay* overlay = m_layerSurface->GetActiveOverlay();
    SurfaceOverlayProperties* p = new SurfaceOverlayProperties( overlay );
    UpdateOverlayProperties( p );
    
    if ( overlay )
    {
      m_histogramWidget->SetInputData( overlay->GetData(), overlay->GetDataSize() );
      
      int nBins = m_histogramWidget->GetNumberOfBins();
      float* fData = new float[ nBins ];
      unsigned char* nColorTable = new unsigned char[ nBins*4 ];
      double range[2];
      overlay->GetRange( range );
      double bin_width = ( range[1] - range[0] ) / nBins;
      int rgb[3];
      double* dColor = m_layerSurface->GetProperties()->GetBinaryColor();
      rgb[0] = (int)( dColor[0] * 255 ); 
      rgb[1] = (int)( dColor[1] * 255 ); 
      rgb[2] = (int)( dColor[2] * 255 ); 
      for ( int i = 0; i < nBins; i++ )
      {
        nColorTable[i*4] = rgb[0];
        nColorTable[i*4+1] = rgb[1];
        nColorTable[i*4+2] = rgb[2];
        nColorTable[i*4+3] = 255;
        
        fData[i] = range[0] + ( i + 0.5 ) * bin_width;
      }
      p->MapOverlayColor( fData, nColorTable, nBins ); 
      m_histogramWidget->SetColorTableData( nColorTable, false );
      delete[] fData;
      delete[] nColorTable;
      
      // rebuild marker lines for display
      MarkerArray markers;
      LineMarker marker;  
      marker.position = 0;
      marker.color = wxColour( 0, 0, 0 );
      marker.style = LineMarker::Dash;
      markers.push_back( marker );
      
      marker.position = p->GetMinPoint();
      marker.color = wxColour( 255, 0, 0 );
      marker.style = LineMarker::Solid;
      markers.push_back( marker );
      
      marker.position = -marker.position;
      markers.push_back( marker );
      
      if ( p->GetColorMethod() == SurfaceOverlayProperties::CM_Piecewise )
      {
        marker.position = p->GetMidPoint();
        marker.color = wxColour( 0, 0, 255 );
        markers.push_back( marker );
        
        marker.position = -marker.position;
        markers.push_back( marker );
      }
      
      marker.position = p->GetMaxPoint();
      marker.color = wxColour( 0, 215, 0 );
      markers.push_back( marker );
      
      marker.position = -marker.position;
      markers.push_back( marker );
      
      m_histogramWidget->SetMarkers( markers );
    }
    
    delete p;
  } 
}

void WindowOverlayConfiguration::OnSliderOpacity( wxScrollEvent& event )
{
  UpdateTextValue( m_textOpacity, event.GetPosition() / 100.0 );
  UpdateGraph();
}

void WindowOverlayConfiguration::OnTextOpacity( wxCommandEvent& event )
{
  if ( m_textOpacity )
  {
    wxString strg = m_textOpacity->GetValue();
    double dValue;
    if ( strg.ToDouble( &dValue ) )
    {
      m_sliderOpacity->SetValue( (int)(dValue * 100) );
      UpdateGraph();
    }
  }
}

void WindowOverlayConfiguration::OnUpdateGraph( wxCommandEvent& event )
{
  UpdateGraph();  
}

void WindowOverlayConfiguration::OnUpdateThreshold( wxCommandEvent& event )
{
  m_textThresholdMid->Enable( m_radioThresholdPiecewise->GetValue() );
  XRCCTRL( *this, "ID_STATIC_MID_POINT", wxStaticText )->Enable( m_radioThresholdPiecewise->GetValue() );
  UpdateGraph();
}


void WindowOverlayConfiguration::OnTextThresholdChanged( wxCommandEvent& event )
{
  UpdateThresholdChanges();
}

void WindowOverlayConfiguration::UpdateThresholdChanges()
{
  if ( !m_radioThresholdPiecewise->GetValue() )   // do not adjust mid point automatically in Piecewise mode
  { 
    double dmin, dmax;
    if ( m_textThresholdMin->GetValue().ToDouble( &dmin ) &&
         m_textThresholdMax->GetValue().ToDouble( &dmax ) )
    {
      UpdateTextValue( m_textThresholdMid, ( dmax + dmin ) / 2 );
    }
  }
  UpdateGraph();
}
