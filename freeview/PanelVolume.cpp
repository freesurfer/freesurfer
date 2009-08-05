/**
 * @file  PanelVolume.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/08/05 17:13:06 $
 *    $Revision: 1.29 $
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

#include "wx/wx.h"

#include "PanelVolume.h"
#include <wx/xrc/xmlres.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerMRI.h"
#include "LayerDTI.h"
#include "LayerPropertiesMRI.h"
#include "LayerPropertiesDTI.h"
#include "wxColorIndicator.h"
#include "DialogEditLookupTable.h"

BEGIN_EVENT_TABLE( PanelVolume, wxPanel )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_VOLUMES" ),        PanelVolume::OnLayerSelectionChanged )
  EVT_CHECKLISTBOX    ( XRCID( "ID_LISTBOX_VOLUMES" ),        PanelVolume::OnLayerVisibilityChanged )
  // EVT_LISTBOX_DCLICK  ( XRCID( "ID_LISTBOX_VOLUMES" ) ),   PanelVolume::OnListDoubleClicked )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_OPACITY" ),         PanelVolume::OnSliderOpacityChanged )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_LOAD" ),            PanelVolume::OnButtonLoad )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_MOVE_UP" ),         PanelVolume::OnButtonMoveUp )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_MOVE_DOWN" ),       PanelVolume::OnButtonMoveDown )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_DELETE" ),          PanelVolume::OnButtonDelete )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_NEW" ),             PanelVolume::OnButtonNew )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_LOAD" ),            PanelVolume::OnButtonLoad )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_SAVE" ),            PanelVolume::OnButtonSave )
  EVT_MENU            ( XRCID( "ID_VOLUME_CLOSE" ),           PanelVolume::OnButtonDelete )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_CLOSE" ),           PanelVolume::OnVolumeCloseUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_MOVE_UP" ),         PanelVolume::OnButtonMoveUp )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_MOVE_UP" ),         PanelVolume::OnMoveUpUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_MOVE_DOWN" ),       PanelVolume::OnButtonMoveDown )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_MOVE_DOWN" ),       PanelVolume::OnMoveDownUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_LOCK" ),            PanelVolume::OnVolumeLock )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_LOCK" ),            PanelVolume::OnVolumeLockUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_COPY_SETTING" ),    PanelVolume::OnVolumeCopySetting )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_COPY_SETTING" ),    PanelVolume::OnVolumeCopySettingUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_PASTE_SETTING" ),   PanelVolume::OnVolumePasteSetting )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_PASTE_SETTING" ),   PanelVolume::OnVolumePasteSettingUpdateUI )
  EVT_MENU            ( XRCID( "ID_VOLUME_PASTE_SETTING_ALL" ), PanelVolume::OnVolumePasteSettingAll )
  EVT_UPDATE_UI       ( XRCID( "ID_VOLUME_PASTE_SETTING_ALL" ), PanelVolume::OnVolumePasteSettingAllUpdateUI )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_CLEAR_BACKGROUND" ), PanelVolume::OnCheckClearBackground )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_SMOOTH" ),        PanelVolume::OnCheckSmooth )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_COLORMAP" ),        PanelVolume::OnChoiceColorMap )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_LUT" ),             PanelVolume::OnChoiceLUT )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_DIRECTION_CODE" ),  PanelVolume::OnChoiceDirectionCode )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_COLORTABLE" ),     PanelVolume::OnColorSelectionChanged )
  EVT_LISTBOX_DCLICK  ( XRCID( "ID_LISTBOX_COLORTABLE" ),     PanelVolume::OnColorSelectionChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_DRAW_VALUE" ),        PanelVolume::OnTextFillValueChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_WINDOW" ),            PanelVolume::OnTextWindowChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_LEVEL" ),             PanelVolume::OnTextLevelChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MIN" ),     PanelVolume::OnTextHeatScaleMinChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MID" ),     PanelVolume::OnTextHeatScaleMidChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MAX" ),     PanelVolume::OnTextHeatScaleMaxChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_OFFSET" ),  PanelVolume::OnTextHeatScaleOffsetChanged )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HEATSCALE_CLEAR_HIGH" ),  PanelVolume::OnCheckHeatScaleClearHigh )
  EVT_TEXT            ( XRCID( "ID_TEXT_JETSCALE_MIN" ),      PanelVolume::OnTextMinJetScaleChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_JETSCALE_MAX" ),      PanelVolume::OnTextMaxJetScaleChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_GRAYSCALE_MIN" ),     PanelVolume::OnTextGrayScaleMin )
  EVT_TEXT            ( XRCID( "ID_TEXT_GRAYSCALE_MAX" ),     PanelVolume::OnTextGrayScaleMax )
  EVT_TEXT            ( XRCID( "ID_TEXT_OPACITY" ),           PanelVolume::OnTextOpacityChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_WINDOW" ),          PanelVolume::OnSliderWindowChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_LEVEL" ),           PanelVolume::OnSliderLevelChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MIN" ),   PanelVolume::OnSliderHeatScaleMinChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MID" ),   PanelVolume::OnSliderHeatScaleMidChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MAX" ),   PanelVolume::OnSliderHeatScaleMaxChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_OFFSET" ),  PanelVolume::OnSliderHeatScaleOffsetChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_JETSCALE_MIN" ),    PanelVolume::OnSliderMinJetScaleChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_JETSCALE_MAX" ),    PanelVolume::OnSliderMaxJetScaleChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_GRAYSCALE_MIN" ),   PanelVolume::OnSliderGrayScaleMin )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_GRAYSCALE_MAX" ),   PanelVolume::OnSliderGrayScaleMax )
  EVT_TEXT            ( XRCID( "ID_TEXT_FRAME" ),             PanelVolume::OnTextFrameChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_FRAME" ),           PanelVolume::OnSliderFrameChanged )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_DISPLAY_VECTOR" ), PanelVolume::OnCheckDisplayVector )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_DISPLAY_TENSOR" ), PanelVolume::OnCheckDisplayTensor )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_INVERSION" ),       PanelVolume::OnChoiceInversion )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_REPRESENTATION" ),  PanelVolume::OnChoiceRepresentation )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_CONTOUR" ),    PanelVolume::OnCheckContour )
  EVT_TEXT            ( XRCID( "ID_TEXT_CONTOUR_MIN" ),    PanelVolume::OnTextContourMin )
  EVT_TEXT            ( XRCID( "ID_TEXT_CONTOUR_MAX" ),    PanelVolume::OnTextContourMax )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_CONTOUR_MIN" ),  PanelVolume::OnSliderContourMin )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_CONTOUR_MAX" ),  PanelVolume::OnSliderContourMax )
END_EVENT_TABLE()


PanelVolume::PanelVolume( wxWindow* parent ) : Listener( "PanelVolume" ), Broadcaster( "PanelVolume" )
{
  m_curCTAB = NULL;
  m_bUINeedUpdate = false;
  m_layerCopied = NULL;

  wxXmlResource::Get()->LoadPanel( this, parent, _("ID_PANEL_VOLUME") );
  m_btnNew              = XRCCTRL( *this, "ID_BUTTON_NEW", wxButton );
  m_btnDelete           = XRCCTRL( *this, "ID_BUTTON_DELETE", wxButton );
  m_btnMoveUp           = XRCCTRL( *this, "ID_BUTTON_MOVE_UP", wxButton );
  m_btnMoveDown         = XRCCTRL( *this, "ID_BUTTON_MOVE_DOWN", wxButton );
  m_btnSave             = XRCCTRL( *this, "ID_BUTTON_SAVE", wxButton );
  m_listBoxLayers       = XRCCTRL( *this, "ID_LISTBOX_VOLUMES", wxCheckListBox );
  m_sliderOpacity       = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
  m_textOpacity         = XRCCTRL( *this, "ID_TEXT_OPACITY", wxTextCtrl );
  m_checkClearBackground =  XRCCTRL( *this, "ID_CHECKBOX_CLEAR_BACKGROUND", wxCheckBox );
  m_checkSmooth         = XRCCTRL( *this, "ID_CHECKBOX_SMOOTH", wxCheckBox );
  m_listColorTable      = XRCCTRL( *this, "ID_LISTBOX_COLORTABLE", wxListBox );
  m_choiceColorMap      = XRCCTRL( *this, "ID_CHOICE_COLORMAP", wxChoice );
  m_choiceLUT           = XRCCTRL( *this, "ID_CHOICE_LUT", wxChoice );
  m_textDrawValue       = XRCCTRL( *this, "ID_TEXT_DRAW_VALUE", wxTextCtrl );
  m_sliderWindow        = XRCCTRL( *this, "ID_SLIDER_WINDOW", wxSlider );
  m_sliderLevel         = XRCCTRL( *this, "ID_SLIDER_LEVEL", wxSlider );
  m_textWindow          = XRCCTRL( *this, "ID_TEXT_WINDOW", wxTextCtrl );
  m_textLevel           = XRCCTRL( *this, "ID_TEXT_LEVEL", wxTextCtrl );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_textHeatScaleMin    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MIN", wxTextCtrl );
  m_textHeatScaleMid    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MID", wxTextCtrl );
  m_textHeatScaleMax    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MAX", wxTextCtrl );
  m_textHeatScaleOffset = XRCCTRL( *this, "ID_TEXT_HEATSCALE_OFFSET", wxTextCtrl );
  m_sliderHeatScaleMin  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MIN", wxSlider );
  m_sliderHeatScaleMid  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MID", wxSlider );
  m_sliderHeatScaleMax  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MAX", wxSlider );
  m_sliderHeatScaleOffset = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_OFFSET", wxSlider );
  m_textJetScaleMin     = XRCCTRL( *this, "ID_TEXT_JETSCALE_MIN", wxTextCtrl );
  m_textJetScaleMax     = XRCCTRL( *this, "ID_TEXT_JETSCALE_MAX", wxTextCtrl );
  m_sliderJetScaleMin   = XRCCTRL( *this, "ID_SLIDER_JETSCALE_MIN", wxSlider );
  m_sliderJetScaleMax   = XRCCTRL( *this, "ID_SLIDER_JETSCALE_MAX", wxSlider );
  m_choiceDirectionCode = XRCCTRL( *this, "ID_CHOICE_DIRECTION_CODE", wxChoice );
  m_colorIndicator      = XRCCTRL( *this, "ID_COLOR_INDICATOR", wxColorIndicator );
  m_textFrame           = XRCCTRL( *this, "ID_TEXT_FRAME", wxTextCtrl );
  m_sliderFrame         = XRCCTRL( *this, "ID_SLIDER_FRAME", wxSlider );
  m_checkDisplayVector  = XRCCTRL( *this, "ID_CHECKBOX_DISPLAY_VECTOR", wxCheckBox );
  m_checkDisplayTensor  = XRCCTRL( *this, "ID_CHECKBOX_DISPLAY_TENSOR", wxCheckBox );
  m_choiceInversion     = XRCCTRL( *this, "ID_CHOICE_INVERSION", wxChoice );
  m_sliderGrayScaleMin  = XRCCTRL( *this, "ID_SLIDER_GRAYSCALE_MIN", wxSlider );
  m_sliderGrayScaleMax  = XRCCTRL( *this, "ID_SLIDER_GRAYSCALE_MAX", wxSlider );
  m_textGrayScaleMin    = XRCCTRL( *this, "ID_TEXT_GRAYSCALE_MIN", wxTextCtrl );
  m_textGrayScaleMax    = XRCCTRL( *this, "ID_TEXT_GRAYSCALE_MAX", wxTextCtrl );
  m_choiceRepresentation  = XRCCTRL( *this, "ID_CHOICE_REPRESENTATION", wxChoice );
  m_choiceMask          = XRCCTRL( *this, "ID_CHOICE_MASK", wxChoice );
  m_checkHeatScaleClearHigh = XRCCTRL( *this, "ID_CHECKBOX_HEATSCALE_CLEAR_HIGH", wxCheckBox );

  // workaround for wxformbuilder bug
  m_colorIndicator->SetSize( 40, wxDefaultCoord );

  m_checkContour =      XRCCTRL( *this, "ID_CHECKBOX_CONTOUR", wxCheckBox );
  m_sliderContourMin =  XRCCTRL( *this, "ID_SLIDER_CONTOUR_MIN", wxSlider );
  m_sliderContourMax =  XRCCTRL( *this, "ID_SLIDER_CONTOUR_MAX", wxSlider );
  m_textContourMin =    XRCCTRL( *this, "ID_TEXT_CONTOUR_MIN", wxTextCtrl );
  m_textContourMax =    XRCCTRL( *this, "ID_TEXT_CONTOUR_MAX", wxTextCtrl );

  m_checkContour->Hide();

  m_luts = MainWindow::GetMainWindowPointer()->GetLUTData();

  m_choiceDirectionCode->Append( _("RAS -> RGB") );
  m_choiceDirectionCode->Append( _("RAS -> RBG") );
  m_choiceDirectionCode->Append( _("RAS -> GRB") );
  m_choiceDirectionCode->Append( _("RAS -> GBR") );
  m_choiceDirectionCode->Append( _("RAS -> BRG") );
  m_choiceDirectionCode->Append( _("RAS -> BGR") );

  m_widgetlistGrayScale.push_back( m_checkClearBackground );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_WINDOW",        wxStaticText ) );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_LEVEL",         wxStaticText ) );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_GRAYSCALE_MIN", wxStaticText ) );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_GRAYSCALE_MAX", wxStaticText ) );
  m_widgetlistGrayScale.push_back( m_textWindow );
  m_widgetlistGrayScale.push_back( m_textLevel );
  m_widgetlistGrayScale.push_back( m_textGrayScaleMin );
  m_widgetlistGrayScale.push_back( m_textGrayScaleMax );
  m_widgetlistGrayScale.push_back( m_sliderWindow );
  m_widgetlistGrayScale.push_back( m_sliderLevel );
  m_widgetlistGrayScale.push_back( m_sliderGrayScaleMin );
  m_widgetlistGrayScale.push_back( m_sliderGrayScaleMax );

  m_widgetlistHeatScale.push_back( m_textHeatScaleMin );
  m_widgetlistHeatScale.push_back( m_textHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_textHeatScaleMax );
  m_widgetlistHeatScale.push_back( m_textHeatScaleOffset );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMin );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMax );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleOffset );
  m_widgetlistHeatScale.push_back( m_checkHeatScaleClearHigh );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MIN",     wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MID",     wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MAX",     wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_OFFSET",  wxStaticText ) );

  m_widgetlistJetScale.push_back( m_textJetScaleMin );
  m_widgetlistJetScale.push_back( m_textJetScaleMax );
  m_widgetlistJetScale.push_back( m_sliderJetScaleMin );
  m_widgetlistJetScale.push_back( m_sliderJetScaleMax );
  m_widgetlistJetScale.push_back( XRCCTRL( *this, "ID_STATIC_JETSCALE_MIN", wxStaticText ) );
  m_widgetlistJetScale.push_back( XRCCTRL( *this, "ID_STATIC_JETSCALE_MAX", wxStaticText ) );

  m_widgetlistLUT.push_back( m_listColorTable );
  m_widgetlistLUT.push_back( m_colorIndicator );
  m_widgetlistLUT.push_back( XRCCTRL( *this, "ID_STATIC_LUT", wxStaticText ) );
  m_widgetlistLUT.push_back( m_choiceLUT );

  m_widgetlistDirectionCode.push_back( m_choiceDirectionCode );
  m_widgetlistDirectionCode.push_back( XRCCTRL( *this, "ID_STATIC_DIRECTION_CODE", wxStaticText ) );

  m_widgetlistFrame.push_back( m_sliderFrame );
  m_widgetlistFrame.push_back( m_textFrame );
  m_widgetlistFrame.push_back( XRCCTRL( *this, "ID_STATIC_FRAME", wxStaticText ) );
  m_widgetlistFrame.push_back( m_checkDisplayVector );
  
  m_widgetlistVector.push_back( XRCCTRL( *this, "ID_STATIC_INVERSION", wxStaticText ) );
  m_widgetlistVector.push_back( m_choiceInversion );
  m_widgetlistVector.push_back( XRCCTRL( *this, "ID_STATIC_REPRESENTATION", wxStaticText ) );
  m_widgetlistVector.push_back( m_choiceRepresentation );
  m_widgetlistVector.push_back( XRCCTRL( *this, "ID_STATIC_MASK", wxStaticText ) );
  m_widgetlistVector.push_back( m_choiceMask );
      
  m_widgetlistContour.push_back( m_sliderContourMin );
  m_widgetlistContour.push_back( m_sliderContourMax );
  m_widgetlistContour.push_back( m_textContourMin );
  m_widgetlistContour.push_back( m_textContourMax );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_MIN", wxStaticText ) );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_MAX", wxStaticText ) );
  
  m_widgetlistEditable.push_back( XRCCTRL( *this, "ID_STATIC_BRUSH_VALUE", wxStaticText ) );
  m_widgetlistEditable.push_back( m_textDrawValue );  
  
  m_widgetlistNormalDisplay.push_back( XRCCTRL( *this, "ID_STATIC_OPACITY", wxStaticText ) );
  m_widgetlistNormalDisplay.push_back( m_sliderOpacity );
  m_widgetlistNormalDisplay.push_back( m_textOpacity );
  m_widgetlistNormalDisplay.push_back( m_checkSmooth );
  m_widgetlistNormalDisplay.push_back( XRCCTRL( *this, "ID_STATIC_COLORMAP", wxStaticText ) );
  m_widgetlistNormalDisplay.push_back( m_choiceColorMap );
  for ( size_t i = 0; i < m_widgetlistGrayScale.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistGrayScale[i] ); 
  for ( size_t i = 0; i < m_widgetlistHeatScale.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistHeatScale[i] );
  for ( size_t i = 0; i < m_widgetlistJetScale.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistJetScale[i] );
  for ( size_t i = 0; i < m_widgetlistLUT.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistLUT[i] );
  for ( size_t i = 0; i < m_widgetlistContour.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistContour[i] );
  for ( size_t i = 0; i < m_widgetlistDirectionCode.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistDirectionCode[i] );
  for ( size_t i = 0; i < m_widgetlistEditable.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistEditable[i] );

  MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->AddListener( this );
  UpdateUI();
}

PanelVolume::~PanelVolume()
{}

void PanelVolume::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// MainWindow* mainwnd = MainWindow::GetMainWindow();
// LayerCollection* lc = mainwnd->GetLayerCollection();
  if ( iMsg == "LayerAdded" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "MRI" ) )
    {
      m_listBoxLayers->Insert( wxString::FromAscii( layer->GetName() ), 0, (void*)layer );
      m_listBoxLayers->Check( 0 );
      m_listBoxLayers->SetSelection( 0 );
      UpdateUI();
    }
  }
  else if ( iMsg == "LayerMoved" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "MRI" ) )
    {
      UpdateLayerList( layer );
      // UpdateUI();
    }
  }
  else if ( iMsg == "LayerModified" || iMsg == /*"WindowLevelChanged"*/ "LayerActorUpdated" ||
            iMsg == "LayerContourShown" )
  {
    UpdateUI();
  }
}

void PanelVolume::OnSliderOpacityChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetOpacity( event.GetPosition() / 100.0 );
  }
}

void PanelVolume::OnTextOpacityChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textOpacity->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetOpacity( dvalue );
    }
  }
}

void PanelVolume::OnLayerSelectionChanged( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  UpdateUI();
}

void PanelVolume::OnLayerVisibilityChanged( wxCommandEvent& event )
{
  int nItem = event.GetInt();
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nItem );
  if ( layer )
  {
    layer->SetVisible( m_listBoxLayers->IsChecked( nItem) );
  }
}

void PanelVolume::OnButtonMoveUp( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

    if ( layer )
      lc->MoveLayerUp( layer );
  }
}

void PanelVolume::OnButtonMoveDown( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

    if ( layer )
      lc->MoveLayerDown( layer );
  }
}

void PanelVolume::UpdateLayerList( Layer* layer )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
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

void PanelVolume::OnButtonDelete( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );
    if ( ((LayerMRI*)layer)->IsModified() )
    {
      wxString msg = _("Volume has been modified. Do you want to close it without saving?");
      wxMessageDialog dlg( this, msg, _("Close"), wxYES_NO | wxICON_QUESTION | wxNO_DEFAULT );
      if ( dlg.ShowModal() != wxID_YES )
        return;
    }
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

void PanelVolume::OnCheckClearBackground( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetClearZero( event.IsChecked() );
  }
}

void PanelVolume::OnCheckSmooth( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetTextureSmoothing( event.IsChecked() ? 1 : 0 );
  }
}

void PanelVolume::OnButtonNew( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->NewVolume();
}

void PanelVolume::OnButtonLoad( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->LoadVolume();
}

void PanelVolume::OnButtonSave( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveVolume();
}

void PanelVolume::OnListDoubleClicked( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      wxTextEntryDialog dlg( this, _("Enter the name of the volume:"), 
                             _("Volume"), 
                             wxString::FromAscii( layer->GetName() ) );
      if ( dlg.ShowModal() == wxID_OK )
      {
        layer->SetName( dlg.GetValue().Trim( true ).Trim( false ).char_str() );
        bool bChecked = m_listBoxLayers->IsChecked( event.GetInt() );
        m_listBoxLayers->SetString( event.GetInt(), wxString::FromAscii( layer->GetName() ) );
        m_listBoxLayers->Check( event.GetInt(), bChecked );
      }
    }
  }
}

void PanelVolume::OnChoiceColorMap( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      int nColorMap = (long)(void*)m_choiceColorMap->GetClientData( event.GetSelection() );
      layer->GetProperties()->SetColorMap( (LayerPropertiesMRI::ColorMapType)nColorMap );
      if ( layer->GetProperties()->GetLUTCTAB() == NULL )
      {
        COLOR_TABLE* ct = m_luts->GetColorTable( 0 );
        layer->GetProperties()->SetLUTCTAB( ct );
      }
    }
  }
  UpdateUI();
}

void PanelVolume::OnChoiceLUT( wxCommandEvent& event )
{
  if ( event.GetSelection() >= m_luts->GetCount() )
  {
    MainWindow::GetMainWindowPointer()->LoadLUT();
  }
  else if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      COLOR_TABLE* ct = m_luts->GetColorTable( event.GetSelection() );
      layer->GetProperties()->SetLUTCTAB( ct );
    }
  }
  UpdateUI();
}

void PanelVolume::PopulateColorTable( COLOR_TABLE* ct )
{
  if ( ct != m_curCTAB )
  {
    m_curCTAB = ct;
    m_listColorTable->Clear();
    m_listColorTable->FitInside();
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    int nValid;
    char name[1000];
    int nSel = -1;
    long nValue = 0;
    m_textDrawValue->GetValue().ToLong( &nValue );
    int nValidCount = 0;
    for ( int i = 0; i < nTotalCount; i++ )
    {
      CTABisEntryValid( ct, i, &nValid );
      if ( nValid )
      {
        CTABcopyName( ct, i, name, 1000 );
        m_listColorTable->Append( wxString::Format( _("%d: %s"), i, name ) );
        if ( nValidCount == nValue )
          nSel = nValidCount;
        nValidCount++;
      }
    }
    if ( nSel >= 0 )
      m_listColorTable->SetSelection( nSel );
  }
}

void PanelVolume::ShowWidgets( std::vector<wxWindow*>& list, bool bShow )
{
  for ( size_t i = 0; i < list.size(); i++ )
  {
    list[i]->Show( bShow );
  }
}

void PanelVolume::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void PanelVolume::DoUpdateUI()
{
	bool bHasVolume = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
	wxWindowList children = GetChildren();
	wxWindowList::iterator it = children.begin(), end = children.end();
	for (; it != end; it++) 
	{
		if ( !(*it)->IsKindOf(CLASSINFO(wxToolBar) ) && *it != m_listBoxLayers ) 
			(*it)->Enable( bHasVolume );
	} 
	
	LayerMRI* layer = NULL;
	LayerPropertiesMRI::ColorMapType nColorMap = LayerPropertiesMRI::NoColorMap;
	if ( bHasVolume )
	{
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
		for ( int i = 0; i < (int)m_listBoxLayers->GetCount() && i < lc->GetNumberOfLayers(); i++ )
			m_listBoxLayers->Check( i, lc->GetLayer( i )->IsVisible() );
		layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
		if ( layer )
		{
			m_listBoxLayers->Check( m_listBoxLayers->GetSelection(), layer->IsVisible() );
			m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
			UpdateTextValue( m_textOpacity, layer->GetProperties()->GetOpacity() );
			m_checkClearBackground->SetValue( layer->GetProperties()->GetClearZero() );
			m_checkSmooth->SetValue( layer->GetProperties()->GetTextureSmoothing() );
			
      // color map settings
			m_choiceColorMap->SetSelection( layer->GetProperties()->GetColorMap() );
      m_choiceLUT->Clear();
      for ( int i = 0; i < m_luts->GetCount(); i++ )
      {
        m_choiceLUT->Append( wxString::FromAscii( m_luts->GetName( i ) ) );
      }
      m_choiceLUT->Append( _("Load lookup table...") );
			m_choiceLUT->SetSelection( m_luts->GetIndex( layer->GetProperties()->GetLUTCTAB() ) );
	
			UpdateTextValue( m_textDrawValue, (double)layer->GetFillValue() );
      double dwindow = layer->GetProperties()->GetWindow();
      double dlevel  = layer->GetProperties()->GetLevel();	
			UpdateTextValue( m_textWindow, dwindow );
      UpdateTextValue( m_textLevel,  dlevel );
      UpdateTextValue( m_textGrayScaleMin, dlevel - dwindow/2 );
      UpdateTextValue( m_textGrayScaleMax, dlevel + dwindow/2 );
			double* windowrange = layer->GetProperties()->GetWindowRange();
			double* levelrange = layer->GetProperties()->GetLevelRange();
			m_sliderWindow->SetValue( (int)( ( layer->GetProperties()->GetWindow() - windowrange[0] ) / ( windowrange[1] - windowrange[0] ) * 100 ) );
			m_sliderLevel->SetValue( (int)( ( layer->GetProperties()->GetLevel() - levelrange[0] ) / ( levelrange[1] - levelrange[0] ) * 100 ) );
      double dminvalue = layer->GetProperties()->GetMinValue();
      double dmaxvalue = layer->GetProperties()->GetMaxValue();
      double range_min = dminvalue - (dmaxvalue-dminvalue)/4;
      double range_max = dmaxvalue + (dmaxvalue-dminvalue)/4;
      m_sliderGrayScaleMin->SetValue( (int)( ( dlevel - dwindow/2 - range_min ) / ( range_max - range_min ) * 100 ) );
      m_sliderGrayScaleMax->SetValue( (int)( ( dlevel + dwindow/2 - range_min ) / ( range_max - range_min ) * 100 ) );
      
			if ( layer->IsTypeOf( "DTI" ) )
        m_textFileName->ChangeValue( wxString::FromAscii( ((LayerDTI*)layer)->GetVectorFileName() ) );
			else
        m_textFileName->ChangeValue( wxString::FromAscii( layer->GetFileName() ) );
			m_textFileName->SetInsertionPointEnd();
			m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
			
			double fMin = layer->GetProperties()->GetMinValue();
			double fMax = layer->GetProperties()->GetMaxValue();
			m_sliderHeatScaleMin->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMinThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderHeatScaleMid->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMidThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderHeatScaleMax->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMaxThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderHeatScaleOffset->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleOffset() + fMax ) / ( fMax + fMax ) * 100 ) );
			UpdateTextValue( m_textHeatScaleMin, layer->GetProperties()->GetHeatScaleMinThreshold() );
			UpdateTextValue( m_textHeatScaleMid, layer->GetProperties()->GetHeatScaleMidThreshold() );
			UpdateTextValue( m_textHeatScaleMax, layer->GetProperties()->GetHeatScaleMaxThreshold() );
			UpdateTextValue( m_textHeatScaleOffset, layer->GetProperties()->GetHeatScaleOffset() );
			
			double fJetMin = fMin - (fMax-fMin)/4;
			double fJetMax = fMax + (fMax-fMin)/4;
			m_sliderJetScaleMin->SetValue( (int)( ( layer->GetProperties()->GetMinJetScaleWindow() - fJetMin ) / ( fJetMax - fJetMin ) * 100 ) );
			m_sliderJetScaleMax->SetValue( (int)( ( layer->GetProperties()->GetMaxJetScaleWindow() - fJetMin ) / ( fJetMax - fJetMin ) * 100 ) );
			UpdateTextValue( m_textJetScaleMin, layer->GetProperties()->GetMinJetScaleWindow() );
			UpdateTextValue( m_textJetScaleMax, layer->GetProperties()->GetMaxJetScaleWindow() );
					
			m_choiceColorMap->Clear();
      m_choiceColorMap->Append( _("Grayscale"), (void*)LayerPropertiesMRI::Grayscale );
      m_choiceColorMap->Append( _("Heat"), (void*)LayerPropertiesMRI::Heat );
      m_choiceColorMap->Append( _("Jet"), (void*)LayerPropertiesMRI::Jet );
			if ( layer->IsTypeOf( "DTI" ) )
			{
        m_choiceColorMap->Append( _("Direction-coded"), (void*)LayerPropertiesMRI::DirectionCoded );
				m_choiceDirectionCode->SetSelection( ((LayerDTI*)layer)->GetProperties()->GetDirectionCode() );
			}
			else
			{
        m_choiceColorMap->Append( _("Lookup Table"), (void*)LayerPropertiesMRI::LUT );		
			}

			nColorMap = layer->GetProperties()->GetColorMap();
			for ( int i = 0; i < (int)m_choiceColorMap->GetCount(); i++ )
			{
				if ( (void*)m_choiceColorMap->GetClientData( i ) == (void*)nColorMap )
				{
					m_choiceColorMap->SetSelection( i );
					break;
				}
			}
			
			int nFrames = layer->GetNumberOfFrames();
			if ( nFrames > 1 )
				m_sliderFrame->SetRange( 1, nFrames );
			m_sliderFrame->SetValue( layer->GetActiveFrame() + 1 );
			UpdateTextValue( m_textFrame, layer->GetActiveFrame() + 1 );
			
			m_checkContour->SetValue( layer->GetProperties()->GetShowAsContour() );
			m_sliderContourMin->SetValue( (int)( ( layer->GetProperties()->GetContourMinThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderContourMax->SetValue( (int)( ( layer->GetProperties()->GetContourMaxThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			UpdateTextValue( m_textContourMin, layer->GetProperties()->GetContourMinThreshold() );
			UpdateTextValue( m_textContourMax, layer->GetProperties()->GetContourMaxThreshold() );
      
      m_choiceRepresentation->Clear();
      if ( layer->GetProperties()->GetDisplayVector() )
      {
        m_choiceRepresentation->Append( _("Simple Line") );
        m_choiceRepresentation->Append( _("3D Bar (slow!)") );
        m_choiceRepresentation->SetSelection( layer->GetProperties()->GetVectorRepresentation() );
        m_choiceInversion->SetSelection( layer->GetProperties()->GetVectorInversion() );
      }
      else if ( layer->GetProperties()->GetDisplayTensor() )
      {
        m_choiceRepresentation->Append( _("Boxoid") );
        m_choiceRepresentation->Append( _("Ellipsoid (Very slow!)") );
        m_choiceRepresentation->SetSelection( layer->GetProperties()->GetTensorRepresentation() );
        m_choiceInversion->SetSelection( layer->GetProperties()->GetTensorInversion() );
      }
		}
	}
	MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
	m_btnNew->Enable( bHasVolume );
	m_btnDelete->Enable( bHasVolume && !mainWnd->IsProcessing() );	
	m_btnMoveUp->Enable( bHasVolume && m_listBoxLayers->GetSelection() != 0 );
	m_btnMoveDown->Enable( bHasVolume && m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
	m_btnSave->Enable( bHasVolume && layer && layer->IsModified() && !mainWnd->IsProcessing() );	
	
  ShowWidgets( m_widgetlistNormalDisplay, layer && 
                                          !layer->GetProperties()->GetDisplayVector() &&
                                          !layer->GetProperties()->GetDisplayTensor());
  if ( !layer || ( !layer->GetProperties()->GetDisplayVector() && !layer->GetProperties()->GetDisplayTensor() ) )
  {
    ShowWidgets( m_widgetlistGrayScale, layer && nColorMap == LayerPropertiesMRI::Grayscale );
    ShowWidgets( m_widgetlistHeatScale, layer && nColorMap == LayerPropertiesMRI::Heat );
    ShowWidgets( m_widgetlistJetScale, layer && nColorMap == LayerPropertiesMRI::Jet );
    ShowWidgets( m_widgetlistLUT, layer && nColorMap == LayerPropertiesMRI::LUT );
    ShowWidgets( m_widgetlistDirectionCode, layer && nColorMap == LayerPropertiesMRI::DirectionCoded );
    ShowWidgets( m_widgetlistEditable, layer && layer->IsEditable() );
  }
  ShowWidgets( m_widgetlistFrame, layer && 
                                  !layer->IsTypeOf( "DTI" ) && 
                                  layer->GetNumberOfFrames() > 1 );
  m_sliderFrame->Enable( layer && 
                         !layer->GetProperties()->GetDisplayVector() && 
                         !layer->GetProperties()->GetDisplayTensor() );
  m_textFrame->Enable( layer && 
                      !layer->GetProperties()->GetDisplayVector() && 
                      !layer->GetProperties()->GetDisplayTensor() );
  m_checkDisplayVector->Show( layer && ( layer->IsTypeOf( "DTI" ) || layer->GetNumberOfFrames() == 3 ) );
  m_checkDisplayVector->SetValue( layer && layer->GetProperties()->GetDisplayVector() );
  m_checkDisplayTensor->Show( layer && layer->GetNumberOfFrames() == 9 );
  m_checkDisplayTensor->SetValue( layer && layer->GetProperties()->GetDisplayTensor() );
  ShowWidgets( m_widgetlistVector, m_checkDisplayVector->IsChecked() || m_checkDisplayTensor->IsChecked() );  
  
//	ShowWidgets( m_widgetlistContour, m_checkContour->IsChecked() );
  ShowWidgets( m_widgetlistContour, false );
  m_checkContour->Show( false /*nColorMap == LayerPropertiesMRI::LUT*/ );
	
	Layout();
	if ( layer && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
	{
		PopulateColorTable( layer->GetProperties()->GetLUTCTAB() );
    
		UpdateColorIndicator();    
    
    m_listColorTable->SetSelection( wxNOT_FOUND );
    for ( int i = 0; i < (int)m_listColorTable->GetCount(); i++ )
    {
      wxString strg = m_listColorTable->GetString( i );
      double dvalue;
      if ( strg.Left( strg.Find( _(":") ) ).ToDouble( &dvalue ) && dvalue == layer->GetFillValue() )
      {
        m_listColorTable->SetSelection( i );
        break;
      }
    }
	}
}

void PanelVolume::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  wxString value_strg = ( (wxString)_("") << dvalue );
  if ( value_strg != ctrl->GetValue() && (value_strg + _(".")) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void PanelVolume::UpdateColorIndicator()
{
  long nIndex = 0;
  if ( m_textDrawValue->GetValue().ToLong( &nIndex ) && m_curCTAB )
  {
    int nr, ng, nb;
    int nValid = 0;
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( m_curCTAB, &nTotalCount );
    if ( nIndex < nTotalCount )
      CTABisEntryValid( m_curCTAB, nIndex, &nValid );
    if ( nValid && CTABrgbAtIndexi( m_curCTAB, nIndex, &nr, &ng, &nb ) == 0 )
    {
      m_colorIndicator->SetColor( wxColour( nr, ng, nb ) );
    }
    else
    {
      m_colorIndicator->SetColor( wxNullColour );
    }
  }
}

void PanelVolume::OnColorSelectionChanged( wxCommandEvent& event )
{
  if ( m_listColorTable->GetSelection() == wxNOT_FOUND )
    return;

  wxString strg = m_listColorTable->GetString( m_listColorTable->GetSelection() );
  double dValue;
  if ( strg.Left( strg.Find( _(":") ) ).ToDouble( &dValue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    layer->SetFillValue( (float)dValue );

    UpdateUI();
  }
}

void PanelVolume::OnTextFillValueChanged( wxCommandEvent& event )
{
  wxString strgvalue = m_textDrawValue->GetValue();
  long nvalue;
  double dvalue;
  if ( strgvalue.ToLong( &nvalue ) )
  {
    m_listColorTable->SetSelection( wxNOT_FOUND );
    for ( int i = 0; i < (int)m_listColorTable->GetCount(); i++ )
    {
      wxString strg = m_listColorTable->GetString( i );
      if ( strg.Left( strg.Find( _(":") ) ).ToDouble( &dvalue ) && dvalue == (double)nvalue )
      {
        m_listColorTable->SetFirstItem( i );
        m_listColorTable->SetSelection( i );
        break;
      }
    }

    UpdateColorIndicator();

    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->SetFillValue( (float)nvalue );
  }
}

void PanelVolume::OnTextWindowChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textWindow->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetWindow() != dvalue )
    {
      // this->BlockListen( true );
      layer->GetProperties()->SetWindow( dvalue );
      // this->BlockListen( false );
    }
  }
}

void PanelVolume::OnTextLevelChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textLevel->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetLevel() != dvalue )
    {
      layer->GetProperties()->SetLevel( dvalue );
    }
  }
}

void PanelVolume::OnSliderWindowChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double* r = layer->GetProperties()->GetWindowRange();
    layer->GetProperties()->SetWindow( (double)m_sliderWindow->GetValue() / 100.0 * ( r[1] - r[0] ) + r[0] );
  }
}

void PanelVolume::OnSliderLevelChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double* r = layer->GetProperties()->GetLevelRange();
    layer->GetProperties()->SetLevel( (double)m_sliderLevel->GetValue() / 100.0 * ( r[1] - r[0] ) + r[0] );
  }
}

void PanelVolume::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxPanel::OnInternalIdle();
}


void PanelVolume::OnTextHeatScaleMinChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMin->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleMinThreshold() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleMinThreshold( dvalue );
    }
  }
}

void PanelVolume::OnTextHeatScaleMidChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMid->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleMidThreshold() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleMidThreshold( dvalue );
    }
  }
}

void PanelVolume::OnTextHeatScaleMaxChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleMaxThreshold() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleMaxThreshold( dvalue );
    }
  }
}

void PanelVolume::OnTextHeatScaleOffsetChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleOffset->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleOffset() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleOffset( dvalue );
    }
  }
}

void PanelVolume::OnSliderHeatScaleMinChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetHeatScaleMinThreshold( (double)m_sliderHeatScaleMin->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderHeatScaleMidChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetHeatScaleMidThreshold( (double)m_sliderHeatScaleMid->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderHeatScaleMaxChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetHeatScaleMaxThreshold( (double)m_sliderHeatScaleMax->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderHeatScaleOffsetChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    // double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetHeatScaleOffset( (double)m_sliderHeatScaleOffset->GetValue() / 100.0 * ( fMax + fMax ) - fMax );
  }
}

void PanelVolume::OnTextMinJetScaleChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textJetScaleMin->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetMinJetScaleWindow() != dvalue )
    {
      layer->GetProperties()->SetMinJetScaleWindow( dvalue );
    }
  }
}

void PanelVolume::OnTextMaxJetScaleChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textJetScaleMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetMaxJetScaleWindow() != dvalue )
    {
      layer->GetProperties()->SetMaxJetScaleWindow( dvalue );
    }
  }
}

void PanelVolume::OnSliderMinJetScaleChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    double fJetMin = fMin - (fMax-fMin)/4;
    double fJetMax = fMax + (fMax-fMin)/4;
    layer->GetProperties()->SetMinJetScaleWindow( (double)m_sliderJetScaleMin->GetValue() / 100.0 * ( fJetMax - fJetMin ) + fJetMin );
  }
}

void PanelVolume::OnSliderMaxJetScaleChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    double fJetMin = fMin - (fMax-fMin)/4;
    double fJetMax = fMax + (fMax-fMin)/4;
    layer->GetProperties()->SetMaxJetScaleWindow( (double)m_sliderJetScaleMax->GetValue() / 100.0 * ( fJetMax - fJetMin ) + fJetMin );
  }
}

void PanelVolume::OnChoiceDirectionCode( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerDTI* layer = ( LayerDTI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->IsTypeOf( "DTI" ) )
    {
      layer->GetProperties()->SetDirectionCode( (LayerPropertiesDTI::DirectionCode)event.GetSelection() );
    }
  }
  UpdateUI();
}

void PanelVolume::OnSliderFrameChanged( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->SetActiveFrame( event.GetInt() - 1 );
  }
}

void PanelVolume::OnTextFrameChanged( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  long value;
  if ( layer && m_textFrame->GetValue().ToLong( &value ) )
  {
    layer->SetActiveFrame( value - 1 );
  }
}

void PanelVolume::PanelVolume::OnVolumeCloseUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND && !MainWindow::GetMainWindowPointer()->IsProcessing() );
}

void PanelVolume::OnMoveUpUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND && m_listBoxLayers->GetSelection() != 0 );
}

void PanelVolume::OnMoveDownUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND &&
                m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
}


void PanelVolume::OnVolumeLock( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->Lock( event.IsChecked() );
  }
}

void PanelVolume::OnVolumeLockUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  event.Check( layer && layer->IsLocked() );
}

void PanelVolume::OnVolumeCopySetting( wxCommandEvent& event )
{
  m_layerCopied = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
}

void PanelVolume::OnVolumeCopySettingUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
}

void PanelVolume::OnVolumePasteSetting( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  layer->GetProperties()->CopySetttings( m_layerCopied->GetProperties() );
}

void PanelVolume::OnVolumePasteSettingUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  event.Enable( lc->Contains( m_layerCopied ) && m_layerCopied != layer );
}


void PanelVolume::OnVolumePasteSettingAll( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    LayerMRI* layer = (LayerMRI*)lc->GetLayer( i );
    if ( layer->GetProperties()->GetColorMap() == m_layerCopied->GetProperties()->GetColorMap() )
      layer->GetProperties()->CopySetttings( m_layerCopied->GetProperties() );
  }
}

void PanelVolume::OnVolumePasteSettingAllUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  event.Enable( lc->Contains( m_layerCopied ) );
}

void PanelVolume::OnCheckContour( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  layer->GetProperties()->SetShowAsContour( event.IsChecked() );
}


void PanelVolume::OnTextContourMin( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textContourMin->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetContourMinThreshold() != dvalue )
    {
      layer->GetProperties()->SetContourMinThreshold( dvalue );
    }
  }
}

void PanelVolume::OnTextContourMax( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textContourMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetContourMaxThreshold() != dvalue )
    {
      layer->GetProperties()->SetContourMaxThreshold( dvalue );
    }
  }
}


void PanelVolume::OnSliderContourMin( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetContourMinThreshold( (double)m_sliderContourMin->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderContourMax( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    layer->GetProperties()->SetContourMaxThreshold( (double)m_sliderContourMax->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderGrayScaleMin( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    double fGrayScaleMin = fMin - (fMax-fMin)/4;
    double fGrayScaleMax = fMax + (fMax-fMin)/4;
    layer->GetProperties()->SetMinGrayscaleWindow( (double)m_sliderGrayScaleMin->GetValue() / 100.0 * ( fGrayScaleMax - fGrayScaleMin ) + fGrayScaleMin );
  }
}

void PanelVolume::OnSliderGrayScaleMax( wxScrollEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    double fMin = layer->GetProperties()->GetMinValue();
    double fMax = layer->GetProperties()->GetMaxValue();
    double fGrayScaleMin = fMin - (fMax-fMin)/4;
    double fGrayScaleMax = fMax + (fMax-fMin)/4;
    layer->GetProperties()->SetMaxGrayscaleWindow( (double)m_sliderGrayScaleMax->GetValue() / 100.0 * ( fGrayScaleMax - fGrayScaleMin ) + fGrayScaleMin );
  }
}

void PanelVolume::OnTextGrayScaleMin( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textGrayScaleMin->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetMinGrayscaleWindow( dvalue );
    }
  }
}


void PanelVolume::OnTextGrayScaleMax( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textGrayScaleMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetMaxGrayscaleWindow( dvalue );
    }
  }
}

void PanelVolume::OnCheckDisplayVector( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->GetProperties()->SetDisplayVector( event.IsChecked() );
    UpdateUI();
  } 
}

void PanelVolume::OnCheckDisplayTensor( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->GetProperties()->SetDisplayTensor( event.IsChecked() );
    UpdateUI();
  } 
}

void PanelVolume::OnChoiceInversion( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    if ( layer->GetProperties()->GetDisplayVector() )
      layer->GetProperties()->SetVectorInversion( event.GetSelection() );
    else if ( layer->GetProperties()->GetDisplayTensor() )
      layer->GetProperties()->SetTensorInversion( event.GetSelection() );
    UpdateUI();
  }   
}

void PanelVolume::OnChoiceRepresentation( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    if ( layer->GetProperties()->GetDisplayVector() )
      layer->GetProperties()->SetVectorRepresentation( event.GetSelection() );
    else if ( layer->GetProperties()->GetDisplayTensor() )
      layer->GetProperties()->SetTensorRepresentation( event.GetSelection() );
    UpdateUI();
  } 
}

void PanelVolume::OnCheckHeatScaleClearHigh( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    layer->GetProperties()->SetHeatScaleClearHigh( event.IsChecked() );
  }
}
