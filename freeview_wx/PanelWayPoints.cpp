/**
 * @file  PanelWayPoints.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:40 $
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

#include "wx/wx.h"
#include <wx/clrpicker.h>
#include <wx/filedlg.h>
#include "PanelWayPoints.h"
#include <wx/xrc/xmlres.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerWayPoints.h"
#include "LayerPropertiesWayPoints.h"
#include "FSSurface.h"
#include "LayerMRI.h"

BEGIN_EVENT_TABLE( PanelWayPoints, wxPanel )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_WAYPOINTS" ),    PanelWayPoints::OnLayerSelectionChanged )
  EVT_CHECKLISTBOX    ( XRCID( "ID_LISTBOX_WAYPOINTS" ),    PanelWayPoints::OnLayerVisibilityChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_OPACITY" ),       PanelWayPoints::OnSliderOpacity )
  EVT_MENU            ( XRCID( "ID_WAYPOINTS_CLOSE" ),      PanelWayPoints::OnWayPointsClose )
  EVT_UPDATE_UI       ( XRCID( "ID_WAYPOINTS_CLOSE" ),      PanelWayPoints::OnWayPointsCloseUpdateUI )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER" ),    PanelWayPoints::OnColorChanged )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER_SPLINE" ),  PanelWayPoints::OnSplineColorChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_RADIUS" ),          PanelWayPoints::OnTextRadiusChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_SPLINE_RADIUS" ),   PanelWayPoints::OnTextSplineRadius )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_COLOR_MAP" ),     PanelWayPoints::OnChoiceColorMap )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_SCALAR_MAP" ),    PanelWayPoints::OnChoiceScalarMap )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MIN" ), PanelWayPoints::OnSliderHeatScaleMin )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MID" ), PanelWayPoints::OnSliderHeatScaleMid )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MAX" ), PanelWayPoints::OnSliderHeatScaleMax )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_OFFSET" ), PanelWayPoints::OnSliderHeatScaleOffset )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MIN" ),   PanelWayPoints::OnTextHeatScaleMin )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MID" ),   PanelWayPoints::OnTextHeatScaleMid )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MAX" ),   PanelWayPoints::OnTextHeatScaleMax )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_OFFSET" ),  PanelWayPoints::OnTextHeatScaleOffset )
  EVT_CHECKBOX        ( XRCID( "ID_CHECK_SHOW_SPLINE" ),    PanelWayPoints::OnCheckShowSpline )
  EVT_CHECKBOX        ( XRCID( "ID_CHECK_SNAP_VOXEL_CENTER" ),    PanelWayPoints::OnCheckSnapToVoxelCenter )

END_EVENT_TABLE()


PanelWayPoints::PanelWayPoints( wxWindow* parent ) :
    Listener( "PanelWayPoints" ),
    Broadcaster( "PanelWayPoints" ),
    m_bUINeedUpdate( false )
{
  wxXmlResource::Get()->LoadPanel( this, parent, _("ID_PANEL_WAYPOINTS") );
  m_listBoxLayers       = XRCCTRL( *this, "ID_LISTBOX_WAYPOINTS", wxCheckListBox );
  m_sliderOpacity       = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );

  m_colorPicker         = XRCCTRL( *this, "ID_COLOR_PICKER", wxColourPickerCtrl );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_colorPickerSpline   = XRCCTRL( *this, "ID_COLOR_PICKER_SPLINE", wxColourPickerCtrl );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_textRadius          = XRCCTRL( *this, "ID_TEXT_RADIUS", wxTextCtrl );
  m_textSplineRadius    = XRCCTRL( *this, "ID_TEXT_SPLINE_RADIUS", wxTextCtrl );
  m_choiceColorMap      = XRCCTRL( *this, "ID_CHOICE_COLOR_MAP", wxChoice );
  m_choiceScalarMap     = XRCCTRL( *this, "ID_CHOICE_SCALAR_MAP", wxChoice );
  m_sliderHeatScaleMin  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MIN", wxSlider );
  m_sliderHeatScaleMid  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MID", wxSlider );
  m_sliderHeatScaleMax  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MAX", wxSlider );
  m_sliderHeatScaleOffset = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_OFFSET", wxSlider );
  m_textHeatScaleMin    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MIN", wxTextCtrl );
  m_textHeatScaleMid    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MID", wxTextCtrl );
  m_textHeatScaleMax    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MAX", wxTextCtrl );
  m_textHeatScaleOffset = XRCCTRL( *this, "ID_TEXT_HEATSCALE_OFFSET", wxTextCtrl );
  m_checkShowSpline     = XRCCTRL( *this, "ID_CHECK_SHOW_SPLINE", wxCheckBox );
  m_checkSnapToVoxelCenter     = XRCCTRL( *this, "ID_CHECK_SNAP_VOXEL_CENTER", wxCheckBox );

  m_widgetlistSolidColor.push_back( m_colorPickerSpline );
  m_widgetlistSolidColor.push_back( XRCCTRL( *this, "ID_STATIC_SPACE_FILLER", wxStaticText ) );

  m_widgetlistHeatScale.push_back( m_choiceScalarMap );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMin );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMax );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleOffset );
  m_widgetlistHeatScale.push_back( m_textHeatScaleMin );
  m_widgetlistHeatScale.push_back( m_textHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_textHeatScaleMax );
  m_widgetlistHeatScale.push_back( m_textHeatScaleOffset );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MIN", wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MID", wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MAX", wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_OFFSET", wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_SCALAR_MAP", wxStaticText ) );
  
  m_widgetlistSpline = m_widgetlistHeatScale;
  m_widgetlistSpline.push_back( m_colorPickerSpline );
  m_widgetlistSpline.push_back( XRCCTRL( *this, "ID_STATIC_SPACE_FILLER", wxStaticText ) );
  m_widgetlistSpline.push_back( XRCCTRL( *this, "ID_STATIC_SPLINE_COLOR", wxStaticText ) );
  m_widgetlistSpline.push_back( XRCCTRL( *this, "ID_CHOICE_COLOR_MAP", wxChoice ) );
  m_widgetlistSpline.push_back( XRCCTRL( *this, "ID_STATIC_SPLINE_RADIUS", wxStaticText ) );
  m_widgetlistSpline.push_back( m_textSplineRadius );

  m_choiceColorMap->Append( _("Solid Color") );
  m_choiceColorMap->Append( _("Heat Scale") );
  m_choiceColorMap->SetSelection( 0 );

  MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->AddListener( this );

  UpdateUI();
}

PanelWayPoints::~PanelWayPoints()
{}

void PanelWayPoints::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// MainWindow* mainwnd = MainWindow::GetMainWindow();
// LayerCollection* lc = mainwnd->GetLayerCollection();
  if ( iMsg == "LayerAdded" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "WayPoints" ) )
    {
      m_listBoxLayers->Insert( wxString::FromAscii( layer->GetName() ), 0, (void*)layer );
      m_listBoxLayers->Check( 0 );
      m_listBoxLayers->SetSelection( 0 );
    }

    UpdateUI();
  }
  else if ( iMsg == "LayerMoved" || iMsg == "LayerNameChanged" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "WayPoints" ) )
    {
      UpdateLayerList( layer );
    }
  }
  else if ( iMsg == "LayerModified" || iMsg == "LayerActorUpdated" )
  {
    UpdateUI();
  }
}

void PanelWayPoints::OnSliderOpacity( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetOpacity( event.GetPosition() / 100.0 );
  }
}

void PanelWayPoints::OnLayerSelectionChanged( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" );
    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  UpdateUI();
}

void PanelWayPoints::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void PanelWayPoints::DoUpdateUI()
{
  bool bHasWayPoints = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  wxWindowList children = GetChildren();
  wxWindowList::iterator it = children.begin(), end = children.end();
  for (; it != end; it++)
  {
    if ( !(*it)->IsKindOf(CLASSINFO(wxToolBar) ) && *it != m_listBoxLayers )
      (*it)->Enable( bHasWayPoints );
  }

  LayerWayPoints* layer = NULL;
  int nColorMap = 0;
  bool bShowSpline = false;
  if ( bHasWayPoints )
  {
    layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
      double* rgb = layer->GetProperties()->GetColor();
      m_colorPicker->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      rgb = layer->GetProperties()->GetSplineColor();
      m_colorPickerSpline->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      m_textFileName->ChangeValue( wxString::FromAscii( layer->GetFileName() ) );
      m_textFileName->SetInsertionPointEnd();
      m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
      m_textRadius->ChangeValue( ( wxString() << layer->GetProperties()->GetRadius() ) );
      m_textSplineRadius->ChangeValue( ( wxString() << layer->GetProperties()->GetSplineRadius() ) );
      nColorMap = layer->GetProperties()->GetColorMap();

      double fMin = layer->GetProperties()->GetScalarMinValue();
      double fMax = layer->GetProperties()->GetScalarMaxValue();
      m_sliderHeatScaleMin->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMin() - fMin ) / ( fMax - fMin ) * 100 ) );
      m_sliderHeatScaleMid->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMid() - fMin ) / ( fMax - fMin ) * 100 ) );
      m_sliderHeatScaleMax->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMax() - fMin ) / ( fMax - fMin ) * 100 ) );
      m_sliderHeatScaleOffset->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleOffset() + fMax ) / ( fMax + fMax ) * 100 ) );
      UpdateTextValue( m_textHeatScaleMin, layer->GetProperties()->GetHeatScaleMin() );
      UpdateTextValue( m_textHeatScaleMid, layer->GetProperties()->GetHeatScaleMid() );
      UpdateTextValue( m_textHeatScaleMax, layer->GetProperties()->GetHeatScaleMax() );
      UpdateTextValue( m_textHeatScaleOffset, layer->GetProperties()->GetHeatScaleOffset() );

      m_choiceColorMap->SetSelection( nColorMap );

      m_choiceScalarMap->Clear();
      std::vector<Layer*> layers = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetLayers();
      int nSel = -1;
      for ( size_t i = 0; i < layers.size(); i++ )
      {
        m_choiceScalarMap->Append( wxString::FromAscii( layers[i]->GetName() ), (void*)layers[i] );
        if ( layer->GetProperties()->GetScalarType() == LayerPropertiesWayPoints::ScalarLayer &&
             layer->GetProperties()->GetScalarLayer() == layers[i] )
        {
          nSel = i;
        }
      }
      std::vector<ScalarValues> svs = layer->GetProperties()->GetScalarSets();
      for ( size_t i = 0; i < svs.size(); i++ )
      {
        m_choiceScalarMap->Append( wxString::FromAscii( svs[i].strName.c_str() ), (void*)NULL );
        if ( layer->GetProperties()->GetScalarType() == LayerPropertiesWayPoints::ScalarSet &&
             layer->GetProperties()->GetScalarSet() == (int)i )
        {
          nSel = i + layers.size();
        }
      }
      m_choiceScalarMap->Append( _("Load..."), (void*)NULL );
      if ( nSel >= 0 )
        m_choiceScalarMap->SetSelection( nSel );
      
      bShowSpline = layer->GetProperties()->GetShowSpline();
      m_checkShowSpline->SetValue( bShowSpline );
      
      m_checkSnapToVoxelCenter->SetValue( layer->GetProperties()->GetSnapToVoxelCenter() );
    }

    // LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" );
    // lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }

// MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  m_colorPicker->Enable( layer );
  m_colorPickerSpline->Enable( layer );

  ShowWidgets( m_widgetlistSpline, bShowSpline );
  ShowWidgets( m_widgetlistSolidColor, bHasWayPoints && nColorMap == LayerPropertiesWayPoints::SolidColor );
  ShowWidgets( m_widgetlistHeatScale, bHasWayPoints && nColorMap == LayerPropertiesWayPoints::HeatScale );
  if ( !bShowSpline )
    ShowWidgets( m_widgetlistSpline, false );

  Layout();
}

void PanelWayPoints::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  wxString value_strg = ( (wxString() << dvalue ) );
  if ( value_strg != ctrl->GetValue() && (value_strg + _(".") ) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void PanelWayPoints::OnLayerVisibilityChanged( wxCommandEvent& event )
{
  int nItem = event.GetInt();
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nItem );
  if ( layer )
  {
    layer->SetVisible( m_listBoxLayers->IsChecked( nItem ) );
  }
}

/*
void PanelWayPoints::OnButtonMoveUp( wxCommandEvent& event )
{
 LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
 int nSel = m_listBoxLayers->GetSelection();
 if ( lc && nSel != wxNOT_FOUND )
 {
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

  if ( layer )
   lc->MoveLayerUp( layer );
 }
}

void PanelWayPoints::OnButtonMoveDown( wxCommandEvent& event )
{
 LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
 int nSel = m_listBoxLayers->GetSelection();
 if ( lc && nSel != wxNOT_FOUND )
 {
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

  if ( layer )
   lc->MoveLayerDown( layer );
 }
}
*/

void PanelWayPoints::UpdateLayerList( Layer* layer )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" );
  int nIndex = lc->GetLayerIndex( layer );
  if ( nIndex != -1 )
  {
    wxString name;
    bool bchecked = false;
    bool bselected = false;
    for ( int i = 0; i < (int)m_listBoxLayers->GetCount(); i++ )
    {
      if ( layer == m_listBoxLayers->GetClientData( i ) )
      {
        name = m_listBoxLayers->GetString( i );
        bchecked = m_listBoxLayers->IsChecked( i );
        bselected = ( m_listBoxLayers->GetSelection() == i );
        m_listBoxLayers->Delete( i );
        break;
      }
    }
    if ( !name.IsEmpty() )
    {
      m_listBoxLayers->Insert( name, nIndex, layer );
      m_listBoxLayers->Check( nIndex, bchecked );
      if ( bselected )
        m_listBoxLayers->SetSelection( nIndex );

      UpdateUI();
    }
  }
}

void PanelWayPoints::OnWayPointsClose( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

    m_listBoxLayers->Delete( nSel );

    if ( (int)m_listBoxLayers->GetCount() > nSel )
      m_listBoxLayers->SetSelection( nSel );
    else if ( nSel - 1 >= 0 )
      m_listBoxLayers->SetSelection( nSel - 1 );

    if ( layer )
      lc->RemoveLayer( layer );

    UpdateUI();
  }
}


void PanelWayPoints::OnWayPointsLoad( wxCommandEvent& event )
{
// MainWindow::GetMainWindowPointer()->LoadSurface();
}

void PanelWayPoints::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxPanel::OnInternalIdle();
}

void PanelWayPoints::OnColorChanged( wxColourPickerEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wxColour c = event.GetColour();
    wp->GetProperties()->SetColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}


void PanelWayPoints::OnSplineColorChanged( wxColourPickerEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wxColour c = event.GetColour();
    wp->GetProperties()->SetSplineColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}


void PanelWayPoints::OnWayPointsCloseUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
}


void PanelWayPoints::OnTextRadiusChanged( wxCommandEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wxString strgvalue = m_textRadius->GetValue();
    double dvalue = 0;
    if ( strgvalue.ToDouble( &dvalue ) && dvalue > 0 )
    {
      this->BlockListen( true );
      wp->GetProperties()->SetRadius( dvalue );
      this->BlockListen( false );
    }
  }
}


void PanelWayPoints::OnTextSplineRadius( wxCommandEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wxString strgvalue = m_textSplineRadius->GetValue();
    double dvalue = 0;
    if ( strgvalue.ToDouble( &dvalue ) && dvalue > 0 )
    {
      this->BlockListen( true );
      wp->GetProperties()->SetSplineRadius( dvalue );
      this->BlockListen( false );
    }
  }
}

void PanelWayPoints::OnChoiceColorMap( wxCommandEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wp->GetProperties()->SetColorMap( event.GetSelection() );
    if ( event.GetSelection() == LayerPropertiesWayPoints::HeatScale &&
         wp->GetProperties()->GetScalarLayer() == NULL )
    {
      wp->GetProperties()->SetScalarLayer( (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetLayer( 0 ) );
    }
  }
}

void PanelWayPoints::OnChoiceScalarMap( wxCommandEvent& event )
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    void* p = event.GetClientData();
    if ( p )
      wp->GetProperties()->SetScalarLayer( (LayerMRI*)p );
    else if ( event.GetSelection() == (int)m_choiceScalarMap->GetCount() - 1 )
      LoadScalarValues();
    else
    {
      int offset = (int)m_choiceScalarMap->GetCount()-1 - wp->GetProperties()->GetNumberOfScalarSets();
      wp->GetProperties()->SetScalarSet( event.GetSelection() - offset );
    }
  }
}

void PanelWayPoints::LoadScalarValues()
{
  LayerWayPoints* wp = ( LayerWayPoints* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->GetActiveLayer();
  if ( wp )
  {
    wxFileDialog dlg( this, _("Open scalar value file"), _(""), _(""),
                      _("Value files (*.*)|*.*"),
                      wxFD_OPEN );
    if ( dlg.ShowModal() == wxID_OK )
    {
      if ( !wp->GetProperties()->LoadScalarsFromFile( dlg.GetPath().char_str() ) )
        cout << "Load scalar values failed" << endl;
    }
  }
}

void PanelWayPoints::OnSliderHeatScaleMin( wxScrollEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetScalarMinValue();
    double fMax = layer->GetProperties()->GetScalarMaxValue();
    layer->GetProperties()->SetHeatScaleMin( (double)m_sliderHeatScaleMin->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelWayPoints::OnSliderHeatScaleMid( wxScrollEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetScalarMinValue();
    double fMax = layer->GetProperties()->GetScalarMaxValue();
    layer->GetProperties()->SetHeatScaleMid( (double)m_sliderHeatScaleMid->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelWayPoints::OnSliderHeatScaleMax( wxScrollEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetScalarMinValue();
    double fMax = layer->GetProperties()->GetScalarMaxValue();
    layer->GetProperties()->SetHeatScaleMax( (double)m_sliderHeatScaleMax->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelWayPoints::OnSliderHeatScaleOffset( wxScrollEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
//  double fMin = layer->GetProperties()->GetScalarMinValue();
    double fMax = layer->GetProperties()->GetScalarMaxValue();
    layer->GetProperties()->SetHeatScaleOffset( (double)m_sliderHeatScaleOffset->GetValue() / 100.0 * ( fMax + fMax ) - fMax );
  }
}

void PanelWayPoints::OnTextHeatScaleMin( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMin->GetValue().ToDouble( &dvalue ) )
  {
    LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleMin( dvalue );
    }
  }
}

void PanelWayPoints::OnTextHeatScaleMid( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMid->GetValue().ToDouble( &dvalue ) )
  {
    LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleMid( dvalue );
    }
  }
}

void PanelWayPoints::OnTextHeatScaleMax( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleMax( dvalue );
    }
  }
}

void PanelWayPoints::OnTextHeatScaleOffset( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleOffset->GetValue().ToDouble( &dvalue ) )
  {
    LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleOffset( dvalue );
    }
  }
}

void PanelWayPoints::OnCheckShowSpline( wxCommandEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->GetProperties()->SetShowSpline( m_checkShowSpline->GetValue() );
  }
}

void PanelWayPoints::ShowWidgets( std::vector<wxWindow*>& list, bool bShow )
{
  for ( size_t i = 0; i < list.size(); i++ )
  {
    list[i]->Show( bShow );
  }
}

void PanelWayPoints::OnCheckSnapToVoxelCenter( wxCommandEvent& event )
{
  LayerWayPoints* layer = ( LayerWayPoints* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->GetProperties()->SetSnapToVoxelCenter( m_checkSnapToVoxelCenter->GetValue() );
  }
}
