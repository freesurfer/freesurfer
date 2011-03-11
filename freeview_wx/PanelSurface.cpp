/**
 * @file  PanelSurface.h
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
#include "PanelSurface.h"
#include <wx/xrc/xmlres.h>
#include <wx/spinctrl.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerSurface.h"
#include "LayerPropertiesSurface.h"
#include "FSSurface.h"
#include "SurfaceOverlay.h"
#include "SurfaceAnnotation.h"
#include "SurfaceLabel.h"
#include "MyUtils.h"

BEGIN_EVENT_TABLE( PanelSurface, wxPanel )
  EVT_MENU            ( XRCID( "ID_SURFACE_CLOSE" ),            PanelSurface::OnSurfaceClose )
  EVT_UPDATE_UI       ( XRCID( "ID_SURFACE_CLOSE" ),            PanelSurface::OnSurfaceCloseUpdateUI )  
  EVT_MENU            ( XRCID( "ID_SURFACE_LOCK" ),             PanelSurface::OnSurfaceLock )
  EVT_UPDATE_UI       ( XRCID( "ID_SURFACE_LOCK" ),             PanelSurface::OnSurfaceLockUpdateUI )  
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_SURFACE" ),          PanelSurface::OnLayerSelectionChanged )
  EVT_CHECKLISTBOX    ( XRCID( "ID_LISTBOX_SURFACE" ),          PanelSurface::OnLayerVisibilityChanged )

  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER" ),        PanelSurface::OnColorChanged )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER_EDGE" ),   PanelSurface::OnEdgeColorChanged )
  EVT_SPINCTRL        ( XRCID( "ID_SPIN_EDGE_THICKNESS" ),      PanelSurface::OnSpinEdgeThickness )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_VECTORS" ),           PanelSurface::OnChoiceVector )
  EVT_SPINCTRL        ( XRCID( "ID_SPIN_VECTOR_POINT_SIZE" ),   PanelSurface::OnSpinVectorPointSize )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER_VECTOR" ), PanelSurface::OnVectorColorChanged )
  
  EVT_CHOICE          ( XRCID( "ID_CHOICE_OVERLAY" ),           PanelSurface::OnChoiceOverlay )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_OVERLAY" ),           PanelSurface::OnButtonConfigureOverlay )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_RENDER_MODE" ),       PanelSurface::OnChoiceRenderMode )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_MESH_COLORMAP" ),     PanelSurface::OnChoiceMeshColorMap )

  EVT_CHOICE          ( XRCID( "ID_CHOICE_ANNOTATION" ),        PanelSurface::OnChoiceAnnotation )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_LABEL" ),             PanelSurface::OnChoiceLabel )
  EVT_COLOURPICKER_CHANGED  ( XRCID( "ID_COLORPICKER_LABEL_COLOR" ),  PanelSurface::OnLabelColorChanged )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_SHOW_VERTICES" ),       PanelSurface::OnCheckShowVertices )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER_VERTEX" ),     PanelSurface::OnColorVertex )
  EVT_SPINCTRL        ( XRCID( "ID_SPIN_VERTEX_POINT_SIZE" ),       PanelSurface::OnSpinVertexPointSize )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HIDE_INFO" ),           PanelSurface::OnCheckHideInfo )

  EVT_CHOICE          ( XRCID( "ID_CHOICE_CURVATURE_MAP" ),         PanelSurface::OnChoiceCurvatureMap )
  EVT_COMMAND_SCROLL_THUMBTRACK ( XRCID( "ID_SLIDER_MID_POINT" ),   PanelSurface::OnSliderMidPointChanging )
  EVT_COMMAND_SCROLL_PAGEDOWN   ( XRCID( "ID_SLIDER_MID_POINT" ),   PanelSurface::OnSliderMidPoint )
  EVT_COMMAND_SCROLL_PAGEUP     ( XRCID( "ID_SLIDER_MID_POINT" ),   PanelSurface::OnSliderMidPoint )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_MID_POINT" ), PanelSurface::OnSliderMidPoint )
  EVT_COMMAND_SCROLL_THUMBTRACK ( XRCID( "ID_SLIDER_SLOPE" ),       PanelSurface::OnSliderSlopeChanging )
  EVT_COMMAND_SCROLL_PAGEDOWN   ( XRCID( "ID_SLIDER_SLOPE" ),       PanelSurface::OnSliderSlope )
  EVT_COMMAND_SCROLL_PAGEUP     ( XRCID( "ID_SLIDER_SLOPE" ),       PanelSurface::OnSliderSlope )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_SLOPE" ),     PanelSurface::OnSliderSlope )
  EVT_TEXT_ENTER      ( XRCID( "ID_TEXT_SLOPE" ),                   PanelSurface::OnTextSlope )
  EVT_TEXT_ENTER      ( XRCID( "ID_TEXT_MID_POINT" ),               PanelSurface::OnTextMidPoint )

  EVT_COMMAND_SCROLL_THUMBTRACK ( XRCID( "ID_SLIDER_OPACITY" ),     PanelSurface::OnSliderOpacityChanging )
  EVT_COMMAND_SCROLL_PAGEDOWN   ( XRCID( "ID_SLIDER_OPACITY" ),     PanelSurface::OnSliderOpacity )
  EVT_COMMAND_SCROLL_PAGEUP     ( XRCID( "ID_SLIDER_OPACITY" ),     PanelSurface::OnSliderOpacity )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_OPACITY" ),   PanelSurface::OnSliderOpacity )
  EVT_TEXT_ENTER                ( XRCID( "ID_TEXT_OPACITY_CTRL" ),       PanelSurface::OnTextOpacity )
  
  EVT_TEXT_ENTER                ( XRCID( "ID_TEXT_POSITION" ),      PanelSurface::OnTextPosition )
END_EVENT_TABLE()


PanelSurface::PanelSurface( wxWindow* parent ) :
    Listener( "PanelSurface" ),
    Broadcaster( "PanelSurface" ),
    m_bUINeedUpdate( false )
{
  wxXmlResource::Get()->LoadPanel( this, parent, _("ID_PANEL_SURFACE") );
  m_listBoxLayers = XRCCTRL( *this, "ID_LISTBOX_SURFACE", wxCheckListBox );
  m_sliderOpacity = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
  m_textOpacity   = XRCCTRL( *this, "ID_TEXT_OPACITY_CTRL", wxTextCtrl );
  
  /*
  m_btnNew   = XRCCTRL( *this, "ID_BUTTON_NEW", wxButton );
  m_btnLoad   = XRCCTRL( *this, "ID_BUTTON_LOAD", wxButton );
  m_btnSave   = XRCCTRL( *this, "ID_BUTTON_SAVE", wxButton );
  m_btnDelete  = XRCCTRL( *this, "ID_BUTTON_DELETE", wxButton );

  m_btnSurfaceMain   = XRCCTRL( *this, "ID_BUTTON_SURFACE_MAIN", wxButton );
  m_btnSurfaceInflated  = XRCCTRL( *this, "ID_BUTTON_SURFACE_INFLATED", wxButton );
  m_btnSurfaceWhite   = XRCCTRL( *this, "ID_BUTTON_SURFACE_WHITE", wxButton );
  m_btnSurfacePial   = XRCCTRL( *this, "ID_BUTTON_SURFACE_PIAL", wxButton );
  m_btnSurfaceOriginal  = XRCCTRL( *this, "ID_BUTTON_SURFACE_ORIGINAL", wxButton );
  */

  m_colorPicker         = XRCCTRL( *this, "ID_COLOR_PICKER",          wxColourPickerCtrl );
  m_colorPickerEdge     = XRCCTRL( *this, "ID_COLOR_PICKER_EDGE",     wxColourPickerCtrl );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME",         wxTextCtrl );
  m_spinEdgeThickness   = XRCCTRL( *this, "ID_SPIN_EDGE_THICKNESS",   wxSpinCtrl );

  m_choiceVector        = XRCCTRL( *this, "ID_CHOICE_VECTORS",          wxChoice );
  m_colorPickerVector   = XRCCTRL( *this, "ID_COLOR_PICKER_VECTOR",     wxColourPickerCtrl );
  m_spinVectorPointSize = XRCCTRL( *this, "ID_SPIN_VECTOR_POINT_SIZE",  wxSpinCtrl );

  m_sliderOpacity       = XRCCTRL( *this, "ID_SLIDER_OPACITY",          wxSlider );

  m_choiceCurvatureMap  = XRCCTRL( *this, "ID_CHOICE_CURVATURE_MAP",    wxChoice );
  m_sliderMidPoint      = XRCCTRL( *this, "ID_SLIDER_MID_POINT",        wxSlider );
  m_sliderSlope         = XRCCTRL( *this, "ID_SLIDER_SLOPE",            wxSlider );
  m_textMidPoint        = XRCCTRL( *this, "ID_TEXT_MID_POINT",          wxTextCtrl );
  m_textSlope           = XRCCTRL( *this, "ID_TEXT_SLOPE",              wxTextCtrl );
  
  m_choiceRenderMode    = XRCCTRL( *this, "ID_CHOICE_RENDER_MODE",      wxChoice );
  m_choiceMeshColorMap  = XRCCTRL( *this, "ID_CHOICE_MESH_COLORMAP",    wxChoice );
  m_checkShowVertices   = XRCCTRL( *this, "ID_CHECKBOX_SHOW_VERTICES",  wxCheckBox );
  m_colorPickerVertex   = XRCCTRL( *this, "ID_COLOR_PICKER_VERTEX",     wxColourPickerCtrl );
  m_spinVertexPointSize = XRCCTRL( *this, "ID_SPIN_VERTEX_POINT_SIZE",  wxSpinCtrl );
  
  m_choiceOverlay       = XRCCTRL( *this, "ID_CHOICE_OVERLAY",        wxChoice );
  m_btnOverlayConfiguration 
      = XRCCTRL( *this, "ID_BUTTON_OVERLAY",        wxButton );
  
  m_choiceAnnotation    = XRCCTRL( *this, "ID_CHOICE_ANNOTATION",     wxChoice );
  m_choiceLabel         = XRCCTRL( *this, "ID_CHOICE_LABEL",          wxChoice );
  m_colorPickerLabel    = XRCCTRL( *this, "ID_COLORPICKER_LABEL_COLOR", wxColourPickerCtrl );
  
  m_textPosition        = XRCCTRL( *this, "ID_TEXT_POSITION",         wxTextCtrl );
  
  m_checkHideInfo       = XRCCTRL( *this, "ID_CHECKBOX_HIDE_INFO",    wxCheckBox );

  m_widgetsSlope.push_back( m_sliderSlope );
  m_widgetsSlope.push_back( m_textSlope );
  m_widgetsSlope.push_back( XRCCTRL( *this, "ID_STATIC_SLOPE", wxStaticText ) );

  m_widgetsMidPoint.push_back( m_sliderMidPoint );
  m_widgetsMidPoint.push_back( m_textMidPoint );
  m_widgetsMidPoint.push_back( XRCCTRL( *this, "ID_STATIC_MID_POINT", wxStaticText ) );

  m_widgetsVector.push_back( m_colorPickerVector );
  m_widgetsVector.push_back( m_spinVectorPointSize );
  m_widgetsVector.push_back( XRCCTRL( *this, "ID_STATIC_VECTOR_COLOR", wxStaticText ) );
  m_widgetsVector.push_back( XRCCTRL( *this, "ID_STATIC_VECTOR_POINT_SIZE", wxStaticText ) );
  
  m_widgetsVertex.push_back( m_colorPickerVertex );
  m_widgetsVertex.push_back( m_spinVertexPointSize );
  m_widgetsVertex.push_back( XRCCTRL( *this, "ID_STATIC_VERTEX_COLOR", wxStaticText ) );
  m_widgetsVertex.push_back( XRCCTRL( *this, "ID_STATIC_VERTEX_POINT_SIZE", wxStaticText ) );

  m_widgetsMesh.push_back( m_choiceMeshColorMap );
  m_widgetsMesh.push_back( XRCCTRL( *this, "ID_STATIC_MESH_COLORMAP", wxStaticText ) );
  
  m_widgetsLabel.push_back( m_colorPickerLabel );
  m_widgetsLabel.push_back( XRCCTRL( *this, "ID_STATIC_LABEL_COLOR", wxStaticText ) );
  
  wxScrolledWindow* sw = XRCCTRL( *this, "ID_SCROLL_WINDOW", wxScrolledWindow );
  sw->SetScrollRate( 5, 5 );
  sw->SetMaxSize( wxSize( -1, 10000 ) );
  
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->AddListener( this );

  UpdateUI( true );
}

PanelSurface::~PanelSurface()
{}

void PanelSurface::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// MainWindow* mainwnd = MainWindow::GetMainWindow();
// LayerCollection* lc = mainwnd->GetLayerCollection();
  if ( iMsg == "LayerAdded" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "Surface" ) )
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
    if ( layer && layer->IsTypeOf( "Surface" ) )
    {
      UpdateLayerList( layer );
    }
  }
  else if ( iMsg == "LayerModified" || iMsg == "LayerActorUpdated" )
  {
    UpdateUI();
  }
}

void PanelSurface::OnLayerSelectionChanged( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  UpdateUI();
}

void PanelSurface::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void PanelSurface::DoUpdateUI()
{
  if ( !IsShown() )
    return;
  
  bool bHasSurface = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  wxWindowList children = XRCCTRL( *this, "ID_SCROLL_WINDOW", wxScrolledWindow )->GetChildren();
  wxWindowList::iterator it = children.begin(), end = children.end();
  for (; it != end; it++)
  {
    if ( !(*it)->IsKindOf(CLASSINFO(wxToolBar) ) && *it != m_listBoxLayers )
      (*it)->Enable( bHasSurface );
  }

  LayerSurface* layer = NULL;
  FSSurface* surf = NULL;
  if ( bHasSurface )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
    for ( int i = 0; i < (int)m_listBoxLayers->GetCount() && i < lc->GetNumberOfLayers(); i++ )
      m_listBoxLayers->Check( i, lc->GetLayer( i )->IsVisible() );
    layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
      UpdateTextValue( m_textOpacity, layer->GetProperties()->GetOpacity() );
      double* rgb = layer->GetProperties()->GetBinaryColor();
      m_colorPicker->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      rgb = layer->GetProperties()->GetEdgeColor();
      m_colorPickerEdge->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      rgb = layer->GetProperties()->GetVectorColor();
      m_colorPickerVector->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      m_textFileName->ChangeValue( wxString::FromAscii( layer->GetFileName() ) );
      m_textFileName->SetInsertionPointEnd();
      m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
      m_spinEdgeThickness->SetValue( layer->GetProperties()->GetEdgeThickness() );
      m_spinVectorPointSize->SetValue( layer->GetProperties()->GetVectorPointSize() );

      m_choiceCurvatureMap->SetSelection( layer->GetProperties()->GetCurvatureMap() );

      UpdateTextValue( m_textMidPoint,  layer->GetProperties()->GetThresholdMidPoint()  );
      UpdateTextValue( m_textSlope,   layer->GetProperties()->GetThresholdSlope()  );
      double range[2];
      layer->GetCurvatureRange( range );
      m_sliderMidPoint->SetValue( (int) ( ( layer->GetProperties()->GetThresholdMidPoint() - range[0] ) / ( range[1] - range[0] ) * 100 ) );
      m_sliderSlope->SetValue( (int) ( layer->GetProperties()->GetThresholdSlope() ) );

      surf = layer->GetSourceSurface();     
       
      m_choiceRenderMode->SetSelection( layer->GetProperties()->GetSurfaceRenderMode() );
      m_choiceMeshColorMap->SetSelection( layer->GetProperties()->GetMeshColorMap() );
      m_checkShowVertices->SetValue( layer->GetProperties()->GetShowVertices() );
      rgb = layer->GetProperties()->GetVertexColor();
      m_colorPickerVertex->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );     
      m_spinVertexPointSize->SetValue( layer->GetProperties()->GetVertexPointSize() );
      
      double* dPos = layer->GetProperties()->GetPosition();
      wxString value_strg = ( (wxString)_("") << dPos [0] << _(" ") << dPos[1] << _(" ") << dPos[2] );
      m_textPosition->ChangeValue( value_strg );
      
      m_checkHideInfo->SetValue( !layer->GetProperties()->GetShowInfo() );
    }

    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }

  m_choiceVector->Clear();
  m_choiceVector->Append( _("Off") );
  if ( surf )
  {
    for ( int i = 0; i < surf->GetNumberOfVectorSets(); i++ )
    {
      m_choiceVector->Append( wxString::FromAscii( surf->GetVectorSetName( i ) ) );
    }
  }
  m_choiceVector->Append( _("Load vector data...") );
  m_choiceVector->SetSelection( surf ? 1 + surf->GetActiveVector() : 0 );
  
  // update overlay controls
  m_choiceOverlay->Clear();
  m_choiceOverlay->Append( _("Off") );
  if ( layer )
  {
    for ( int i = 0; i < layer->GetNumberOfOverlays(); i++ )
    {
      m_choiceOverlay->Append( wxString::FromAscii( layer->GetOverlay( i )->GetName() ) );
    }
  }
  m_choiceOverlay->Append( _("Load from file...") );
  m_choiceOverlay->SetSelection( layer ? 1 + layer->GetActiveOverlayIndex() : 0 );
  
  m_btnOverlayConfiguration->Show( layer && layer->GetActiveOverlayIndex() >= 0 );
  
  // update annotation controls
  m_choiceAnnotation->Clear();
  m_choiceAnnotation->Append( _("Off") );
  if ( layer )
  {
    for ( int i = 0; i < layer->GetNumberOfAnnotations(); i++ )
    {
      m_choiceAnnotation->Append( wxString::FromAscii( layer->GetAnnotation( i )->GetName() ) );
    }
  }
  m_choiceAnnotation->Append( _("Load from file...") );
  m_choiceAnnotation->SetSelection( layer ? 1 + layer->GetActiveAnnotationIndex() : 0 );
  
  // update label controls
  m_choiceLabel->Clear();
  if ( layer )
  {
    if ( layer->GetNumberOfLabels() > 0 )
    {
      for ( int i = 0; i < layer->GetNumberOfLabels(); i++ )
      {
        m_choiceLabel->Append( wxString::FromAscii( layer->GetLabel( i )->GetName() ) );
      }
    }
    else
      m_choiceLabel->Append( _("None") );
  }
  m_choiceLabel->Append( _("Load from file...") );
  if ( layer && layer->GetActiveLabelIndex() >= 0 )
    m_choiceLabel->SetSelection( layer->GetActiveLabelIndex() );
  else
    m_choiceLabel->SetSelection( 0 );  
  if ( layer && layer->GetActiveLabel() )
  {
    double* rgb = layer->GetActiveLabel()->GetColor();
    m_colorPickerLabel->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );  
  }   

  int nCurvatureMap = layer ? layer->GetProperties()->GetCurvatureMap() : 0;
  for ( size_t i = 0; i < m_widgetsMidPoint.size(); i++ )
  {
    m_widgetsMidPoint[i]->Show( nCurvatureMap != LayerPropertiesSurface::CM_Off );
  }
  for ( size_t i = 0; i < m_widgetsSlope.size(); i++ )
  {
    m_widgetsSlope[i]->Show( nCurvatureMap == LayerPropertiesSurface::CM_Threshold );
  }
  for ( size_t i = 0; i < m_widgetsVector.size(); i++ )
  {
    m_widgetsVector[i]->Show( m_choiceVector->GetSelection() > 0 );
  }
  for ( size_t i = 0; i < m_widgetsVertex.size(); i++ )
  {
    m_widgetsVertex[i]->Show( m_checkShowVertices->GetValue() );
  }
  for ( size_t i = 0; i < m_widgetsMesh.size(); i++ )
  {
    m_widgetsMesh[i]->Show( layer && layer->GetProperties()->GetSurfaceRenderMode() != LayerPropertiesSurface::SM_Surface );
  }
  for ( size_t i = 0; i < m_widgetsLabel.size(); i++ )
  {
    m_widgetsLabel[i]->Show( layer && layer->GetActiveLabelIndex() >= 0 );
  }
  m_colorPicker->Enable( layer ); // && nCurvatureMap != LayerPropertiesSurface::CM_Threshold );

  // hack to force resize of these controls in scrolled window
  for ( size_t i = 0; i < m_widgetsResize.size(); i++ )
  {
    wxSize sz = m_widgetsResize[i]->GetMinSize();
    m_widgetsResize[i]->SetMinSize( wxSize( 100, sz.GetHeight() ) );
  }
 
  Layout();
}

void PanelSurface::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  wxString value_strg = ( (wxString)_("") << dvalue );
  if ( value_strg != ctrl->GetValue() && (value_strg + _(".")) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void PanelSurface::OnLayerVisibilityChanged( wxCommandEvent& event )
{
  int nItem = event.GetInt();
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nItem );
  if ( layer )
  {
    layer->SetVisible( m_listBoxLayers->IsChecked( nItem ) );
  }
}

/*
void PanelSurface::OnButtonMoveUp( wxCommandEvent& event )
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

void PanelSurface::OnButtonMoveDown( wxCommandEvent& event )
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

void PanelSurface::OnSurfaceLock( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->Lock( event.IsChecked() );
    }
  }
}

void PanelSurface::OnSurfaceLockUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    event.Check( layer && layer->IsLocked() );
  }
}

void PanelSurface::OnSurfaceClose( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );    

    m_listBoxLayers->Delete( nSel );   

    if ( (int)m_listBoxLayers->GetCount() > nSel )
      m_listBoxLayers->SetSelection( nSel );
    else if ( nSel >= 1 )
      m_listBoxLayers->SetSelection( nSel - 1 );        

    if ( layer )
      lc->RemoveLayer( layer );

    UpdateUI();
  }
}

void PanelSurface::OnSurfaceCloseUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
}

void PanelSurface::UpdateLayerList( Layer* layer )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
  int nIndex = lc->GetLayerIndex( layer );
  std::vector<Layer*> layers = lc->GetLayers();
  if ( nIndex != -1 )
  {
    m_listBoxLayers->Clear();
    int nSel = 0;
    for ( size_t i = 0; i < layers.size(); i++ )
    {
      m_listBoxLayers->Append( wxString::FromAscii( layers[i]->GetName() ), layers[i] );
      m_listBoxLayers->Check( i, layers[i]->IsVisible() );
      if ( lc->GetActiveLayer() == layers[i] )
        nSel = i;
    }
    m_listBoxLayers->SetSelection( nSel );

    UpdateUI();
  }
}

void PanelSurface::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxPanel::OnInternalIdle();
}

void PanelSurface::OnColorChanged( wxColourPickerEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    wxColour c = event.GetColour();
    surf->GetProperties()->SetBinaryColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

void PanelSurface::OnEdgeColorChanged( wxColourPickerEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    wxColour c = event.GetColour();
    surf->GetProperties()->SetEdgeColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

void PanelSurface::OnSpinEdgeThickness( wxSpinEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetEdgeThickness( event.GetInt() );
  }
}

void PanelSurface::OnCheckShowVertices( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->ShowVertices( event.IsChecked() );
  }
}


void PanelSurface::OnColorVertex( wxColourPickerEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    wxColour c = event.GetColour();
    surf->GetProperties()->SetVertexColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}


void PanelSurface::OnSpinVertexPointSize( wxSpinEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetVertexPointSize( event.GetInt() );
  }
}


void PanelSurface::OnChoiceVector( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    int nVector = event.GetSelection() - 1;
    if ( nVector < surf->GetNumberOfVectorSets() )
      surf->SetActiveVector( nVector );
    else
    {
      // load new
      MainWindow::GetMainWindowPointer()->LoadSurfaceVector();
      UpdateUI();
    }
  }
}

void PanelSurface::OnSpinVectorPointSize( wxSpinEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetVectorPointSize( event.GetInt() );
  }
}

void PanelSurface::OnVectorColorChanged( wxColourPickerEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    wxColour c = event.GetColour();
    surf->GetProperties()->SetVectorColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

void PanelSurface::OnChoiceCurvatureMap( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    if ( event.GetSelection() < 3 )
      surf->GetProperties()->SetCurvatureMap( event.GetSelection() );
    else
    {
      // load new curvature maps
      MainWindow::GetMainWindowPointer()->LoadSurfaceCurvature();
    }
    UpdateUI();
  }
}

void PanelSurface::OnTextMidPoint( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textMidPoint->GetValue().ToDouble( &dvalue ) )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetThresholdMidPoint() != dvalue )
    {
      layer->GetProperties()->SetThresholdMidPoint( dvalue );
    }
  }
}

void PanelSurface::OnSliderSlope( wxScrollEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    double fMin = 0;
    double fMax = 100;
    surf->GetProperties()->SetThresholdSlope( (double)m_sliderSlope->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelSurface::OnTextSlope( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textSlope->GetValue().ToDouble( &dvalue ) )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetThresholdSlope() != dvalue )
    {
      layer->GetProperties()->SetThresholdSlope( dvalue );
    }
  }
}

void PanelSurface::OnSliderMidPointChanging( wxScrollEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {  
    double range[2];
    surf->GetCurvatureRange( range );
    UpdateTextValue( m_textMidPoint, (double)m_sliderMidPoint->GetValue() / 100.0 * ( range[1] - range[0] ) + range[0] );
  }
}

void PanelSurface::OnSliderMidPoint( wxScrollEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    double range[2];
    surf->GetCurvatureRange( range );
    surf->GetProperties()->SetThresholdMidPoint( (double)m_sliderMidPoint->GetValue() / 100.0 * ( range[1] - range[0] ) + range[0] );
  }
}

void PanelSurface::OnSliderSlopeChanging( wxScrollEvent& event )
{
  double fMin = 0;
  double fMax = 100;
  UpdateTextValue( m_textSlope, (double)m_sliderSlope->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
}


void PanelSurface::OnSliderOpacity( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetOpacity( event.GetPosition() / 100.0 );
  }
}

void PanelSurface::OnSliderOpacityChanging( wxScrollEvent& event )
{
  double fMin = 0;
  double fMax = 1;
  UpdateTextValue( m_textOpacity, (double)m_sliderOpacity->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
}

void PanelSurface::OnTextOpacity( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textOpacity->GetValue().ToDouble( &dvalue ) )
  {
    LayerSurface* layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetOpacity() != dvalue )
    {
      layer->GetProperties()->SetOpacity( dvalue );
    }
  }
}

void PanelSurface::OnButtonConfigureOverlay( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->ConfigureOverlay();
}

void PanelSurface::OnChoiceOverlay( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    int nSel = event.GetSelection() - 1;
    if ( nSel < surf->GetNumberOfOverlays() )
    {
      surf->SetActiveOverlay( nSel );
    }
    else 
    {
      // load new overlay map
       MainWindow::GetMainWindowPointer()->LoadSurfaceOverlay();
    }

    UpdateUI();
  }  
}

void PanelSurface::OnChoiceRenderMode( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetSurfaceRenderMode( event.GetSelection() );
  }  
}


void PanelSurface::OnChoiceMeshColorMap( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetMeshColorMap( event.GetSelection() );
  }  
}


void PanelSurface::OnChoiceAnnotation( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    int nSel = event.GetSelection() - 1;
    if ( nSel < surf->GetNumberOfAnnotations() )
    {
      surf->SetActiveAnnotation( nSel );
    }
    else
    {
      // load new overlay map
      MainWindow::GetMainWindowPointer()->LoadSurfaceAnnotation();
    }
    UpdateUI();
  }  
}

void PanelSurface::OnChoiceLabel( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    int nSel = event.GetSelection();
    if ( nSel >= surf->GetNumberOfLabels() )
    {
      MainWindow::GetMainWindowPointer()->LoadSurfaceLabel();
    }
    else
    {
      surf->SetActiveLabel( nSel );
    }
    UpdateUI();
  }  
}

void PanelSurface::OnLabelColorChanged( wxColourPickerEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    wxColour c = event.GetColour();
    surf->GetActiveLabel()->SetColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

void PanelSurface::OnTextPosition( wxCommandEvent& event )
{
  wxArrayString sa = MyUtils::SplitString( m_textPosition->GetValue(), _(",") );
  if ( sa.Count() < 3 )
    sa = MyUtils::SplitString( m_textPosition->GetValue(), _(" ") );
  
  double pos[3];
  if ( sa.Count() < 3 || !sa[0].ToDouble( pos ) || !sa[1].ToDouble( pos+1 ) || !sa[2].ToDouble( pos+2 ) )
  {
    cerr << "Please enter 3 numbers." << endl;
    return;
  }
      
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetPosition( pos );
    UpdateUI();
  }  
}

void PanelSurface::OnCheckHideInfo( wxCommandEvent& event )
{
  LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf )
  {
    surf->GetProperties()->SetShowInfo( !event.IsChecked() );
  }  
}

