/**
 * @file  PanelVolume.h
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

#include "PanelVolume.h"
#include <wx/clrpicker.h>
#include <wx/xrc/xmlres.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerMRI.h"
#include "LayerDTI.h"
#include "LayerSurface.h"
#include "LayerPropertiesMRI.h"
#include "LayerPropertiesDTI.h"
#include "wxColorIndicator.h"
#include "FSVolume.h"
#include "DialogEditLookupTable.h"

BEGIN_EVENT_TABLE( PanelVolume, wxPanel )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_VOLUMES" ),        PanelVolume::OnLayerSelectionChanged )
  EVT_CHECKLISTBOX    ( XRCID( "ID_LISTBOX_VOLUMES" ),        PanelVolume::OnLayerVisibilityChanged )
  EVT_LISTBOX_DCLICK  ( XRCID( "ID_LISTBOX_VOLUMES" ),        PanelVolume::OnListDoubleClicked )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_OPACITY" ),         PanelVolume::OnSliderOpacityChanged )
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
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_UPSAMPLE" ),      PanelVolume::OnCheckUpsample )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_COLORMAP" ),        PanelVolume::OnChoiceColorMap )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_LUT" ),             PanelVolume::OnChoiceLUT )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_DIRECTION_CODE" ),  PanelVolume::OnChoiceDirectionCode )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_COLORTABLE" ),     PanelVolume::OnColorSelectionChanged )
  EVT_LISTBOX_DCLICK  ( XRCID( "ID_LISTBOX_COLORTABLE" ),     PanelVolume::OnColorSelectionChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_DRAW_VALUE" ),        PanelVolume::OnTextFillValueChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_WINDOW" ),            PanelVolume::OnTextWindowChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_LEVEL" ),             PanelVolume::OnTextLevelChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_COLORMAP_MIN" ),      PanelVolume::OnTextColorMapMinChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_COLORMAP_MAX" ),      PanelVolume::OnTextColorMapMaxChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_MID" ),     PanelVolume::OnTextHeatScaleMidChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_HEATSCALE_OFFSET" ),  PanelVolume::OnTextHeatScaleOffsetChanged )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HEATSCALE_CLEAR_HIGH" ),  PanelVolume::OnCheckHeatScaleClearHigh )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HEATSCALE_TRUNCATE" ),    PanelVolume::OnCheckHeatScaleTruncate )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HEATSCALE_INVERT" ),      PanelVolume::OnCheckHeatScaleInvert )
  EVT_TEXT            ( XRCID( "ID_TEXT_OPACITY_CTRL" ),           PanelVolume::OnTextOpacityChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_WINDOW" ),          PanelVolume::OnSliderWindowChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_LEVEL" ),           PanelVolume::OnSliderLevelChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_COLORMAP_MIN" ),    PanelVolume::OnSliderColorMapMinChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_COLORMAP_MAX" ),    PanelVolume::OnSliderColorMapMaxChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_MID" ),   PanelVolume::OnSliderHeatScaleMidChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_HEATSCALE_OFFSET" ),  PanelVolume::OnSliderHeatScaleOffsetChanged )
  EVT_TEXT            ( XRCID( "ID_TEXT_FRAME_CTRL" ),             PanelVolume::OnTextFrameChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_FRAME" ),           PanelVolume::OnSliderFrameChanged )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_DISPLAY_VECTOR" ), PanelVolume::OnCheckDisplayVector )
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_DISPLAY_TENSOR" ), PanelVolume::OnCheckDisplayTensor )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_INVERSION" ),       PanelVolume::OnChoiceInversion )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_REPRESENTATION" ),  PanelVolume::OnChoiceRepresentation )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_LABEL_OUTLINE" ),      PanelVolume::OnCheckShowLabelOutline )
  
  EVT_CHECKBOX              ( XRCID( "ID_CHECKBOX_CONTOUR" ),     PanelVolume::OnCheckContour )
  EVT_TEXT_ENTER            ( XRCID( "ID_TEXT_CONTOUR_MIN" ),     PanelVolume::OnTextContourMin )
  EVT_TEXT_ENTER            ( XRCID( "ID_TEXT_CONTOUR_MAX" ),     PanelVolume::OnTextContourMax )
  EVT_COLOURPICKER_CHANGED  ( XRCID( "ID_COLORPICKER_CONTOUR" ),  PanelVolume::OnColorContour )
  EVT_CHECKBOX        ( XRCID( "ID_CHECK_USE_IMAGE_COLORMAP" ),   PanelVolume::OnCheckUseImageColorMap )
  EVT_CHECKBOX        ( XRCID( "ID_CHECK_CONTOUR_EXTRACT_ALL" ),   PanelVolume::OnCheckContourExtractAll )
  
  EVT_TEXT_ENTER      ( XRCID( "ID_TEXT_SMOOTH_ITERATIONS" ),     PanelVolume::OnTextContourSmooth )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_SAVE_ISOSURFACE" ),     PanelVolume::OnButtonSaveContour )
  
  EVT_CHECKBOX        ( XRCID( "ID_CHECKBOX_HIDE_INFO" ),         PanelVolume::OnCheckHideInfo )
  EVT_CHOICE          ( XRCID( "ID_CHOICE_UPSAMPLE_METHOD" ),     PanelVolume::OnChoiceUpSampleMethod )
  
  EVT_COMMAND_SCROLL_PAGEDOWN     ( XRCID( "ID_SLIDER_CONTOUR_MIN" ),  PanelVolume::OnSliderContourMin )
  EVT_COMMAND_SCROLL_PAGEUP       ( XRCID( "ID_SLIDER_CONTOUR_MIN" ),  PanelVolume::OnSliderContourMin )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_CONTOUR_MIN" ),  PanelVolume::OnSliderContourMin )
  EVT_COMMAND_SCROLL_THUMBTRACK   ( XRCID( "ID_SLIDER_CONTOUR_MIN" ),  PanelVolume::OnSliderContourMinChanging )
  EVT_COMMAND_SCROLL_PAGEDOWN     ( XRCID( "ID_SLIDER_CONTOUR_MAX" ),  PanelVolume::OnSliderContourMax )
  EVT_COMMAND_SCROLL_PAGEUP       ( XRCID( "ID_SLIDER_CONTOUR_MAX" ),  PanelVolume::OnSliderContourMax )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_CONTOUR_MAX" ),  PanelVolume::OnSliderContourMax )
  EVT_COMMAND_SCROLL_THUMBTRACK   ( XRCID( "ID_SLIDER_CONTOUR_MAX" ),  PanelVolume::OnSliderContourMaxChanging )
  EVT_COMMAND_SCROLL_PAGEDOWN     ( XRCID( "ID_SLIDER_SMOOTH_ITERATIONS" ),  PanelVolume::OnSliderContourSmooth )
  EVT_COMMAND_SCROLL_PAGEUP       ( XRCID( "ID_SLIDER_SMOOTH_ITERATIONS" ),  PanelVolume::OnSliderContourSmooth )
  EVT_COMMAND_SCROLL_THUMBRELEASE ( XRCID( "ID_SLIDER_SMOOTH_ITERATIONS" ),  PanelVolume::OnSliderContourSmooth )
  EVT_COMMAND_SCROLL_THUMBTRACK   ( XRCID( "ID_SLIDER_SMOOTH_ITERATIONS" ),  PanelVolume::OnSliderContourSmoothChanging )
END_EVENT_TABLE()


PanelVolume::PanelVolume(wxWindow* parent) : 
  Listener( "PanelVolume" ), 
  Broadcaster( "PanelVolume" ),
  m_listBoxLayers( NULL ),
  m_textOpacity( NULL )
{
  InitWidgetsFromXRC( parent );
}

void PanelVolume::InitWidgetsFromXRC( wxWindow* parent )
{
  wxXmlResource::Get()->LoadPanel( this, parent, wxT("ID_PANEL_VOLUME"));
  m_listBoxLayers       = XRCCTRL( *this, "ID_LISTBOX_VOLUMES", wxCheckListBox );
  m_checkClearBackground =  XRCCTRL( *this, "ID_CHECKBOX_CLEAR_BACKGROUND", wxCheckBox );
  m_checkSmooth         = XRCCTRL( *this, "ID_CHECKBOX_SMOOTH", wxCheckBox );
  m_checkUpsample       = XRCCTRL( *this, "ID_CHECKBOX_UPSAMPLE", wxCheckBox );
  m_listColorTable      = XRCCTRL( *this, "ID_LISTBOX_COLORTABLE", wxListBox );
  m_choiceColorMap      = XRCCTRL( *this, "ID_CHOICE_COLORMAP", wxChoice );
  m_choiceLUT           = XRCCTRL( *this, "ID_CHOICE_LUT", wxChoice );
  m_textOpacity         = XRCCTRL( *this, "ID_TEXT_OPACITY_CTRL", wxTextCtrl );
  m_textDrawValue       = XRCCTRL( *this, "ID_TEXT_DRAW_VALUE", wxTextCtrl );
  m_sliderOpacity       = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
  m_sliderWindow        = XRCCTRL( *this, "ID_SLIDER_WINDOW", wxSlider );
  m_sliderLevel         = XRCCTRL( *this, "ID_SLIDER_LEVEL", wxSlider );
  m_textWindow          = XRCCTRL( *this, "ID_TEXT_WINDOW", wxTextCtrl );
  m_textLevel           = XRCCTRL( *this, "ID_TEXT_LEVEL", wxTextCtrl );
  m_textFileName        = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_textColorMapMin     = XRCCTRL( *this, "ID_TEXT_COLORMAP_MIN", wxTextCtrl );
  m_textColorMapMax     = XRCCTRL( *this, "ID_TEXT_COLORMAP_MAX", wxTextCtrl );
  m_textHeatScaleMid    = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MID", wxTextCtrl );
  m_textHeatScaleOffset = XRCCTRL( *this, "ID_TEXT_HEATSCALE_OFFSET", wxTextCtrl );
  m_sliderColorMapMin   = XRCCTRL( *this, "ID_SLIDER_COLORMAP_MIN", wxSlider );
  m_sliderColorMapMax   = XRCCTRL( *this, "ID_SLIDER_COLORMAP_MAX", wxSlider );
  m_sliderHeatScaleMid  = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MID", wxSlider );
  m_sliderHeatScaleOffset = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_OFFSET", wxSlider );
  m_choiceDirectionCode = XRCCTRL( *this, "ID_CHOICE_DIRECTION_CODE", wxChoice );
//  m_colorIndicator      = XRCCTRL( *this, "ID_COLOR_INDICATOR", wxColorIndicator );
  m_textFrame           = XRCCTRL( *this, "ID_TEXT_FRAME_CTRL", wxTextCtrl );
  m_sliderFrame         = XRCCTRL( *this, "ID_SLIDER_FRAME", wxSlider );
  m_checkDisplayVector  = XRCCTRL( *this, "ID_CHECKBOX_DISPLAY_VECTOR", wxCheckBox );
  m_checkDisplayTensor  = XRCCTRL( *this, "ID_CHECKBOX_DISPLAY_TENSOR", wxCheckBox );
  m_choiceInversion     = XRCCTRL( *this, "ID_CHOICE_INVERSION", wxChoice );
  m_choiceRepresentation  = XRCCTRL( *this, "ID_CHOICE_REPRESENTATION", wxChoice );
  m_choiceMask          = XRCCTRL( *this, "ID_CHOICE_MASK", wxChoice );
  m_checkHeatScaleClearHigh = XRCCTRL( *this, "ID_CHECKBOX_HEATSCALE_CLEAR_HIGH", wxCheckBox );
  m_checkHeatScaleTruncate  = XRCCTRL( *this, "ID_CHECKBOX_HEATSCALE_TRUNCATE", wxCheckBox );
  m_checkHeatScaleInvert    = XRCCTRL( *this, "ID_CHECKBOX_HEATSCALE_INVERT", wxCheckBox );
  
  m_choiceUpSampleMethod  = XRCCTRL( *this, "ID_CHOICE_UPSAMPLE_METHOD", wxChoice );

  // workaround for wxWidgets 2.9.1 & wxformbuilder bug
  wxPanel* cpanel = XRCCTRL( *this, "ID_COLOR_INDICATOR_PANEL", wxPanel );
  if ( cpanel )
  { 
    wxBoxSizer* sizer = new wxBoxSizer( wxHORIZONTAL );
    cpanel->SetSizer( sizer );
    m_colorIndicator = new wxColorIndicator( cpanel );
    sizer->Add( m_colorIndicator, 1, wxEXPAND|wxALIGN_CENTER_VERTICAL|wxALL, 5 );
  }
  
  m_checkShowLabelOutline   =      XRCCTRL( *this, "ID_CHECKBOX_LABEL_OUTLINE", wxCheckBox );
  
  m_checkContour =      XRCCTRL( *this, "ID_CHECKBOX_CONTOUR", wxCheckBox );
  m_sliderContourMin =  XRCCTRL( *this, "ID_SLIDER_CONTOUR_MIN", wxSlider );
  m_sliderContourMax =  XRCCTRL( *this, "ID_SLIDER_CONTOUR_MAX", wxSlider );
  m_textContourMin =    XRCCTRL( *this, "ID_TEXT_CONTOUR_MIN", wxTextCtrl );
  m_textContourMax =    XRCCTRL( *this, "ID_TEXT_CONTOUR_MAX", wxTextCtrl );
  m_checkUseImageColorMap   = XRCCTRL( *this, "ID_CHECK_USE_IMAGE_COLORMAP", wxCheckBox );
  m_checkContourExtractAll  = XRCCTRL( *this, "ID_CHECK_CONTOUR_EXTRACT_ALL", wxCheckBox );
  m_colorpickerContour      = XRCCTRL( *this, "ID_COLORPICKER_CONTOUR", wxColourPickerCtrl );
  m_sliderContourSmooth =  XRCCTRL( *this, "ID_SLIDER_SMOOTH_ITERATIONS", wxSlider );
  m_textContourSmooth =  XRCCTRL( *this, "ID_TEXT_SMOOTH_ITERATIONS", wxTextCtrl );
  m_btnSaveContour    = XRCCTRL( *this, "ID_BUTTON_SAVE_ISOSURFACE", wxButton );
//  m_checkContour->Hide();

  m_checkHideInfo = XRCCTRL( *this, "ID_CHECKBOX_HIDE_INFO", wxCheckBox );
  
  m_luts = MainWindow::GetMainWindowPointer()->GetLUTData();

  m_widgetlistGrayScale.push_back( m_checkClearBackground );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_WINDOW",        wxStaticText ) );
  m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_LEVEL",         wxStaticText ) );
  m_widgetlistGrayScale.push_back( m_textWindow );
  m_widgetlistGrayScale.push_back( m_textLevel );
  m_widgetlistGrayScale.push_back( m_sliderWindow );
  m_widgetlistGrayScale.push_back( m_sliderLevel );

  m_widgetlistHeatScale.push_back( m_textHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_textHeatScaleOffset );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleMid );
  m_widgetlistHeatScale.push_back( m_sliderHeatScaleOffset );
  m_widgetlistHeatScale.push_back( m_checkHeatScaleClearHigh );
  m_widgetlistHeatScale.push_back( m_checkHeatScaleTruncate );
  m_widgetlistHeatScale.push_back( m_checkHeatScaleInvert );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MID",     wxStaticText ) );
  m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_OFFSET",  wxStaticText ) );

  m_widgetlistGenericColorMap.push_back( m_textColorMapMin );
  m_widgetlistGenericColorMap.push_back( m_textColorMapMax );
  m_widgetlistGenericColorMap.push_back( m_sliderColorMapMin );
  m_widgetlistGenericColorMap.push_back( m_sliderColorMapMax );
  m_widgetlistGenericColorMap.push_back( XRCCTRL( *this, "ID_STATIC_COLORMAP_MIN", wxStaticText ) );
  m_widgetlistGenericColorMap.push_back( XRCCTRL( *this, "ID_STATIC_COLORMAP_MAX", wxStaticText ) );

  m_widgetlistLUT.push_back( m_listColorTable );
  m_widgetlistLUT.push_back( XRCCTRL( *this, "ID_STATIC_LUT", wxStaticText ) );
  m_widgetlistLUT.push_back( m_choiceLUT );
  m_widgetlistLUT.push_back( m_colorIndicator );

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
  m_widgetlistContour.push_back( m_checkUseImageColorMap );
  m_widgetlistContour.push_back( m_checkContourExtractAll );
  m_widgetlistContour.push_back( m_colorpickerContour );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_MIN", wxStaticText ) );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_MAX", wxStaticText ) );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_COLOR", wxStaticText ) );
  m_widgetlistContour.push_back( m_sliderContourSmooth );
  m_widgetlistContour.push_back( m_textContourSmooth );
  m_widgetlistContour.push_back( m_btnSaveContour );
  m_widgetlistContour.push_back( XRCCTRL( *this, "ID_STATIC_SMOOTH_SURFACE", wxStaticText ) );
  
  m_widgetlistEditable.push_back( XRCCTRL( *this, "ID_STATIC_BRUSH_VALUE", wxStaticText ) );
  m_widgetlistEditable.push_back( m_textDrawValue );  
  
  m_widgetlistNormalDisplay.push_back( XRCCTRL( *this, "ID_STATIC_OPACITY", wxStaticText ) );
  m_widgetlistNormalDisplay.push_back( m_sliderOpacity );
  m_widgetlistNormalDisplay.push_back( m_textOpacity );
  m_widgetlistNormalDisplay.push_back( m_checkSmooth );
  m_widgetlistNormalDisplay.push_back( m_checkUpsample );
  m_widgetlistNormalDisplay.push_back( XRCCTRL( *this, "ID_STATIC_COLORMAP", wxStaticText ) );
  m_widgetlistNormalDisplay.push_back( m_choiceColorMap );
  for ( size_t i = 0; i < m_widgetlistGrayScale.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistGrayScale[i] ); 
  for ( size_t i = 0; i < m_widgetlistHeatScale.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistHeatScale[i] );
  for ( size_t i = 0; i < m_widgetlistGenericColorMap.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistGenericColorMap[i] );
  for ( size_t i = 0; i < m_widgetlistLUT.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistLUT[i] );
  for ( size_t i = 0; i < m_widgetlistContour.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistContour[i] );
  for ( size_t i = 0; i < m_widgetlistDirectionCode.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistDirectionCode[i] );
  for ( size_t i = 0; i < m_widgetlistEditable.size(); i++ )
    m_widgetlistNormalDisplay.push_back( m_widgetlistEditable[i] );

  wxScrolledWindow* sw = XRCCTRL( *this, "ID_SCROLL_WINDOW", wxScrolledWindow );
  sw->SetScrollRate( 5, 5 );
  sw->SetMaxSize( wxSize( -1, 10000 ) );
  
  m_widgetlistResize.push_back( m_choiceColorMap );
  m_widgetlistResize.push_back( m_choiceLUT );
  m_widgetlistResize.push_back( m_choiceDirectionCode );
  m_widgetlistResize.push_back( m_choiceInversion );
  m_widgetlistResize.push_back( m_choiceRepresentation );
  m_widgetlistResize.push_back( m_choiceMask );
  m_widgetlistResize.push_back( m_listBoxLayers );
  m_widgetlistResize.push_back( m_listColorTable );
  
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->AddListener( this );
  
  m_curCTAB = NULL;
  m_bUINeedUpdate = false;
  m_layerCopied = NULL;
  
  UpdateUI( true );
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
  else if ( iMsg == "LayerMoved" || iMsg == "LayerNameChanged" )
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
  if ( m_textOpacity && m_listBoxLayers && 
       m_textOpacity->GetValue().ToDouble( &dvalue ) && 
       m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
  LayerCollection* lc_surf = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
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
    for ( int i = 0; i < lc_surf->GetNumberOfLayers(); i++ )
    {
      if ( ( (LayerSurface*)lc_surf->GetLayer( i ) )->GetRefVolume() == layer )
      {
        wxString msg = _("One of the surfaces is using this volume as coordinate reference. You can not close it before closing the surface first.");
        wxMessageDialog dlg( this, msg, _("Close"), wxOK );
        dlg.ShowModal();
        return;
      }
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

void PanelVolume::OnCheckUpsample( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetUpSampleMethod( event.IsChecked() ? LayerPropertiesMRI::UM_BiLinear : LayerPropertiesMRI::UM_None );
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
      MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->MoveToTop( layer );
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
      if ( layer->GetProperties()->GetLUTCTAB() == NULL )
      {
        COLOR_TABLE* ct = m_luts->GetColorTable( 0 );
        layer->GetProperties()->SetLUTCTAB( ct );
      }
      layer->GetProperties()->SetColorMap( nColorMap );
    }
  }
  UpdateUI();
}

void PanelVolume::OnChoiceLUT( wxCommandEvent& event )
{
  if ( event.GetSelection() == (int)m_choiceLUT->GetCount()-1 )
  {
    MainWindow::GetMainWindowPointer()->LoadLUT();
  }
  else if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      if ( event.GetSelection() < m_luts->GetCount() )
      {
        COLOR_TABLE* ct = m_luts->GetColorTable( event.GetSelection() );
        layer->GetProperties()->SetLUTCTAB( ct );
      }
      else
        layer->GetProperties()->SetLUTCTAB( layer->GetEmbeddedColorTable() );
    }
  }
  UpdateUI();
}

void PanelVolume::PopulateColorTable( COLOR_TABLE* ct )
{
  if ( ct && ( ct != m_curCTAB || m_bColorTableNeedReset ) )
  {
    m_curCTAB = ct;
    m_listColorTable->Clear();
    m_listColorTable->FitInside();
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    int nValid = 0;
    char name[1000];
    int nSel = -1;
    long nValue = 0;
    if ( !m_textDrawValue->GetValue().ToLong( &nValue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
    {
      LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
      if ( layer )
        nValue = (int)layer->GetFillValue();
    }
      
    int nValidCount = 0;
    for ( int i = 0; i < nTotalCount; i++ )
    {
      CTABisEntryValid( ct, i, &nValid );
      if ( nValid )
      {
        CTABcopyName( ct, i, name, 1000 );
        m_listColorTable->Append( wxString::Format( _("%d: %s"), i, name ) );
        if ( i == nValue )
          nSel = nValidCount;
        nValidCount++;
      }
    }
    if ( nSel >= 0 )
      m_listColorTable->SetSelection( nSel );
    
    m_bColorTableNeedReset = false;
  }
}

void PanelVolume::SearchInColorTable( const wxString& search_strg )
{
  if ( !m_curCTAB )
    return;
  
  m_listColorTable->Clear();
  m_listColorTable->FitInside();
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_curCTAB, &nTotalCount );
  char name[1000];
  int nValid;
  int nCount = 0;
  int nSel = -1;
  int nValue = 0;
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
    nValue = (int)layer->GetFillValue();
  for ( int i = 0; i < nTotalCount; i++ )
  {
    CTABisEntryValid( m_curCTAB, i, &nValid );
    if ( nValid )
    {
      CTABcopyName( m_curCTAB, i, name, 1000 );
      wxString strName = name;
      strName.MakeLower();
      if ( strName.Find( search_strg.Lower() ) != wxNOT_FOUND )
      {
        m_listColorTable->Append( wxString::Format( _("%d: %s"), i, name ) );
        if ( i == nValue )
          nSel = nCount;
        nCount++;
      }
    }
  }
  if ( nSel >= 0 )
    m_listColorTable->SetSelection( nSel );
  
  m_bColorTableNeedReset = true;
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
  if ( !IsShown() )
    return;
  
	bool bHasVolume = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  wxWindowList children = XRCCTRL( *this, "ID_SCROLL_WINDOW", wxScrolledWindow )->GetChildren();
	wxWindowList::iterator it = children.begin(), end = children.end();
	for (; it != end; it++) 
	{
		if ( !(*it)->IsKindOf(CLASSINFO(wxToolBar) ) && *it != m_listBoxLayers ) 
			(*it)->Enable( bHasVolume );
	} 
	
	LayerMRI* layer = NULL;
	int nColorMap = LayerPropertiesMRI::NoColorMap;
  if ( bHasVolume && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
	{
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
		for ( int i = 0; i < (int)m_listBoxLayers->GetCount() && i < lc->GetNumberOfLayers(); i++ )
			m_listBoxLayers->Check( i, lc->GetLayer( i )->IsVisible() );
		layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
		if ( layer )
		{
      nColorMap = layer->GetProperties()->GetColorMap();
			m_listBoxLayers->Check( m_listBoxLayers->GetSelection(), layer->IsVisible() );
			m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
			UpdateTextValue( m_textOpacity, layer->GetProperties()->GetOpacity() );
			m_checkClearBackground->SetValue( layer->GetProperties()->GetClearZero() );
      m_checkSmooth->SetValue( layer->GetProperties()->GetTextureSmoothing() );
      m_checkUpsample->SetValue( layer->GetProperties()->GetUpSampleMethod() != 0 );
      m_checkUpsample->Enable( layer->GetProperties()->GetColorMap() != LayerPropertiesMRI::LUT );
			
      // color map settings
			m_choiceColorMap->SetSelection( layer->GetProperties()->GetColorMap() );
      m_choiceLUT->Clear();
      for ( int i = 0; i < m_luts->GetCount(); i++ )
      {
        m_choiceLUT->Append( wxString::FromAscii( m_luts->GetName( i ) ) );
      }
      if ( layer->GetEmbeddedColorTable() )
        m_choiceLUT->Append( _("Embedded") );
      m_choiceLUT->Append( _("Load lookup table...") );
      int nSel = m_luts->GetIndex( layer->GetProperties()->GetLUTCTAB() );
      m_choiceLUT->SetSelection( nSel >= 0 ? nSel : m_luts->GetCount() );
	
			UpdateTextValue( m_textDrawValue, (double)layer->GetFillValue() );
      double dwindow = layer->GetProperties()->GetWindow();
      double dlevel  = layer->GetProperties()->GetLevel();	
			UpdateTextValue( m_textWindow, dwindow );
      UpdateTextValue( m_textLevel,  dlevel );
      double dMinTh = dlevel - dwindow/2;
      double dMaxTh = dlevel + dwindow/2;
      double dminvalue = layer->GetProperties()->GetMinValue();
      double dmaxvalue = layer->GetProperties()->GetMaxValue();
      double range_min = dminvalue - (dmaxvalue-dminvalue)/4;
      double range_max = dmaxvalue + (dmaxvalue-dminvalue)/4;
      if ( nColorMap == LayerPropertiesMRI::Heat )
      {
        dMinTh = layer->GetProperties()->GetHeatScaleMinThreshold();
        dMaxTh = layer->GetProperties()->GetHeatScaleMaxThreshold();
        range_min = dminvalue;
        range_max = dmaxvalue;
      }
      else if ( nColorMap != LayerPropertiesMRI::Grayscale && nColorMap != LayerPropertiesMRI::LUT)
      {
        dMinTh = layer->GetProperties()->GetMinGenericThreshold();
        dMaxTh = layer->GetProperties()->GetMaxGenericThreshold();
      }
			double* windowrange = layer->GetProperties()->GetWindowRange();
			double* levelrange = layer->GetProperties()->GetLevelRange();
			m_sliderWindow->SetValue( (int)( ( dwindow - windowrange[0] ) / ( windowrange[1] - windowrange[0] ) * 100 ) );
			m_sliderLevel->SetValue( (int)( ( dlevel - levelrange[0] ) / ( levelrange[1] - levelrange[0] ) * 100 ) );
             
      UpdateTextValue( m_textColorMapMin, dMinTh );
      UpdateTextValue( m_textColorMapMax, dMaxTh );
      m_sliderColorMapMin->SetValue( (int)( ( dMinTh - range_min ) / ( range_max - range_min ) * 100 ) );
      m_sliderColorMapMax->SetValue( (int)( ( dMaxTh - range_min ) / ( range_max - range_min ) * 100 ) );
      
			if ( layer->IsTypeOf( "DTI" ) )
        m_textFileName->ChangeValue( wxString::FromAscii( ((LayerDTI*)layer)->GetVectorFileName() ) );
			else
        m_textFileName->ChangeValue( wxString::FromAscii( layer->GetFileName() ) );
			m_textFileName->SetInsertionPointEnd();
			m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
			
      m_sliderHeatScaleMid->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMidThreshold() - range_min ) / 
          ( range_max - range_min ) * 100 ) );
      m_sliderHeatScaleOffset->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleOffset() + range_max ) / 
          ( range_max + range_max ) * 100 ) );
			UpdateTextValue( m_textHeatScaleMid, layer->GetProperties()->GetHeatScaleMidThreshold() );
      UpdateTextValue( m_textHeatScaleOffset, layer->GetProperties()->GetHeatScaleOffset() );
      m_checkHeatScaleClearHigh->SetValue( layer->GetProperties()->GetHeatScaleClearHigh() );
      m_checkHeatScaleTruncate->SetValue( layer->GetProperties()->GetHeatScaleTruncate() );
      m_checkHeatScaleInvert->SetValue( layer->GetProperties()->GetHeatScaleInvert() );
						
			m_choiceColorMap->Clear();
      m_choiceColorMap->Append( _("Grayscale"), (void*)LayerPropertiesMRI::Grayscale );
			if ( layer->IsTypeOf( "DTI" ) )
			{
        m_choiceColorMap->Append( _("Direction-coded"), (void*)LayerPropertiesMRI::DirectionCoded );
				m_choiceDirectionCode->SetSelection( ((LayerDTI*)layer)->GetProperties()->GetDirectionCode() );
			}
			else
			{
        m_choiceColorMap->Append( _("Lookup Table"), (void*)LayerPropertiesMRI::LUT );		
			}

      m_choiceColorMap->Append( _("Heat"), (void*)LayerPropertiesMRI::Heat );
      m_choiceColorMap->Append( _("Jet"), (void*)LayerPropertiesMRI::Jet );
      m_choiceColorMap->Append( _("GE Color"), (void*)LayerPropertiesMRI::GEColor );
      m_choiceColorMap->Append( _("NIH"), (void*)LayerPropertiesMRI::NIH );
			for ( int i = 0; i < (int)m_choiceColorMap->GetCount(); i++ )
			{
        if ( ((long)(void*)m_choiceColorMap->GetClientData( i )) == ((long)nColorMap) )
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
      m_sliderContourMin->SetValue( (int)( ( layer->GetProperties()->GetContourMinThreshold() - dminvalue ) / ( dmaxvalue - dminvalue ) * 100 ) );
      m_sliderContourMax->SetValue( (int)( ( layer->GetProperties()->GetContourMaxThreshold() - dminvalue ) / ( dmaxvalue - dminvalue ) * 100 ) );
			UpdateTextValue( m_textContourMin, layer->GetProperties()->GetContourMinThreshold() );
			UpdateTextValue( m_textContourMax, layer->GetProperties()->GetContourMaxThreshold() );
      m_checkUseImageColorMap->SetValue( layer->GetProperties()->GetContourUseImageColorMap() );
      m_checkContourExtractAll->SetValue( layer->GetProperties()->GetContourExtractAllRegions() );
      m_sliderContourSmooth->SetValue( layer->GetProperties()->GetContourSmoothIterations() );
      UpdateTextValue( m_textContourSmooth, layer->GetProperties()->GetContourSmoothIterations() );
      
      m_colorpickerContour->Enable( !layer->GetProperties()->GetContourUseImageColorMap() );
      double rgb[3];
      layer->GetProperties()->GetContourColor( rgb );
      m_colorpickerContour->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      
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
      
      m_checkHideInfo->SetValue( !layer->GetProperties()->GetShowInfo() );
      
      m_checkShowLabelOutline->SetValue( layer->GetProperties()->GetShowLabelOutline() );
      m_checkShowLabelOutline->Show( nColorMap == LayerPropertiesMRI::LUT );
      
      m_choiceUpSampleMethod->SetSelection( layer->GetProperties()->GetUpSampleMethod() );
		}
	}
//	MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  ShowWidgets( m_widgetlistNormalDisplay, layer && 
                                          !layer->GetProperties()->GetDisplayVector() &&
                                          !layer->GetProperties()->GetDisplayTensor());
  if ( !layer || ( !layer->GetProperties()->GetDisplayVector() && !layer->GetProperties()->GetDisplayTensor() ) )
  {
    ShowWidgets( m_widgetlistGrayScale, layer && nColorMap == LayerPropertiesMRI::Grayscale );
    ShowWidgets( m_widgetlistHeatScale, layer && nColorMap == LayerPropertiesMRI::Heat );
    ShowWidgets( m_widgetlistGenericColorMap, layer && nColorMap != LayerPropertiesMRI::LUT && 
                                              nColorMap != LayerPropertiesMRI::DirectionCoded );
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
  
	ShowWidgets( m_widgetlistContour, m_checkContour->IsChecked() );
//  ShowWidgets( m_widgetlistContour, false );
//  m_checkContour->Show( false /*nColorMap == LayerPropertiesMRI::LUT*/ );

	if ( layer && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
	{
    if ( m_curCTAB != layer->GetProperties()->GetLUTCTAB() )
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
  
  // hack to force resize of these controls in scrolled window
  for ( size_t i = 0; i < m_widgetlistResize.size(); i++ )
  {
    wxSize sz = m_widgetlistResize[i]->GetMinSize();
    m_widgetlistResize[i]->SetMinSize( wxSize( 100, sz.GetHeight() ) );
  }
 
  Layout();
}

void PanelVolume::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  double val;
  if ( ctrl->GetValue().ToDouble(&val) && val == dvalue )
    return;
  
  wxString value_strg;
  if (dvalue < 1)
    value_strg = ( (wxString)_("") << dvalue );
  else
  {
    value_strg.Printf(_("%f"), dvalue);
    while(value_strg[value_strg.Length()-1] == '0')
      value_strg = value_strg.Left(value_strg.Length()-1);
    if (value_strg[value_strg.Length()-1] == '.')
      value_strg = value_strg.Left(value_strg.Length()-1);
  }
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
  wxString strgvalue = m_textDrawValue->GetValue().Trim( true ).Trim( false );
  long nvalue;
  double dvalue;
  if ( strgvalue.IsEmpty() )
  {
    PopulateColorTable( m_curCTAB );
  }
  else if ( strgvalue.ToLong( &nvalue ) )
  {
    m_listColorTable->SetSelection( wxNOT_FOUND );
    bool bFound = false;
    for ( int i = 0; i < (int)m_listColorTable->GetCount(); i++ )
    {
      wxString strg = m_listColorTable->GetString( i );
      if ( strg.Left( strg.Find( _(":") ) ).ToDouble( &dvalue ) && dvalue == (double)nvalue )
      {
        m_listColorTable->SetFirstItem( i );
        m_listColorTable->SetSelection( i );
        bFound = true;
        break;
      }
    }

    UpdateColorIndicator();

    LayerMRI* layer = NULL;
    if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
      layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->SetFillValue( (float)nvalue );
    
    if ( !bFound )
    {
      PopulateColorTable( m_curCTAB );
    }
  }
  else if ( !strgvalue.IsEmpty() ) // do a search
  {
    SearchInColorTable( strgvalue );
  }
}

void PanelVolume::OnTextWindowChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textWindow->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
  if ( m_textLevel->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double* r = layer->GetProperties()->GetWindowRange();
      layer->GetProperties()->SetWindow( (double)m_sliderWindow->GetValue() / 100.0 * ( r[1] - r[0] ) + r[0] );
    }
  }
}

void PanelVolume::OnSliderLevelChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double* r = layer->GetProperties()->GetLevelRange();
      layer->GetProperties()->SetLevel( (double)m_sliderLevel->GetValue() / 100.0 * ( r[1] - r[0] ) + r[0] );
    }
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


void PanelVolume::OnTextColorMapMinChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textColorMapMin->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      switch ( layer->GetProperties()->GetColorMap() )
      {
        case LayerPropertiesMRI::Grayscale:
          layer->GetProperties()->SetMinGrayscaleWindow( dvalue );
          break;
        case LayerPropertiesMRI::Heat:
          layer->GetProperties()->SetHeatScaleMinThreshold( dvalue );
          break;
        default:
          layer->GetProperties()->SetMinGenericThreshold( dvalue );
          break;
      }
    }
  }
}

void PanelVolume::OnTextColorMapMaxChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textColorMapMax->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      switch ( layer->GetProperties()->GetColorMap() )
      {
        case LayerPropertiesMRI::Grayscale:
          layer->GetProperties()->SetMaxGrayscaleWindow( dvalue );
          break;
        case LayerPropertiesMRI::Heat:
          layer->GetProperties()->SetHeatScaleMaxThreshold( dvalue );
          break;
        default:
          layer->GetProperties()->SetMaxGenericThreshold( dvalue );
          break;
      }
    }
  }
}

void PanelVolume::OnTextHeatScaleMidChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleMid->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleMidThreshold() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleMidThreshold( dvalue );
    }
  }
}

void PanelVolume::OnTextHeatScaleOffsetChanged( wxCommandEvent& event )
{
  double dvalue;
  if ( m_textHeatScaleOffset->GetValue().ToDouble( &dvalue ) && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetHeatScaleOffset() != dvalue )
    {
      layer->GetProperties()->SetHeatScaleOffset( dvalue );
    }
  }
}

void PanelVolume::OnSliderColorMapMinChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      double fScaleMin = fMin - (fMax-fMin)/4;
      double fScaleMax = fMax + (fMax-fMin)/4;
      switch ( layer->GetProperties()->GetColorMap() )
      {
        case LayerPropertiesMRI::Grayscale:
          layer->GetProperties()->SetMinGrayscaleWindow( (double)m_sliderColorMapMin->GetValue() / 
              100.0 * ( fScaleMax - fScaleMin ) + fScaleMin );
          break;
        case LayerPropertiesMRI::Heat:       
          layer->GetProperties()->SetHeatScaleMinThreshold( (double)m_sliderColorMapMin->GetValue() / 
              100.0 * ( fMax - fMin ) + fMin );
          break;
        default:
          layer->GetProperties()->SetMinGenericThreshold( (double)m_sliderColorMapMin->GetValue() / 
              100.0 * ( fScaleMax - fScaleMin ) + fScaleMin  );
          break;
      }
    }
  }
}


void PanelVolume::OnSliderColorMapMaxChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      double fScaleMin = fMin - (fMax-fMin)/4;
      double fScaleMax = fMax + (fMax-fMin)/4;
      switch ( layer->GetProperties()->GetColorMap() )
      {
        case LayerPropertiesMRI::Grayscale:
          layer->GetProperties()->SetMaxGrayscaleWindow( (double)m_sliderColorMapMax->GetValue() / 
              100.0 * ( fScaleMax - fScaleMin ) + fScaleMin );
          break;
        case LayerPropertiesMRI::Heat:       
          layer->GetProperties()->SetHeatScaleMaxThreshold( (double)m_sliderColorMapMax->GetValue() / 
              100.0 * ( fMax - fMin ) + fMin );
          break;
        default:
          layer->GetProperties()->SetMaxGenericThreshold( (double)m_sliderColorMapMax->GetValue() / 
              100.0 * ( fScaleMax - fScaleMin ) + fScaleMin  );
          break;
      }
    }
  }
}

void PanelVolume::OnSliderHeatScaleMidChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      layer->GetProperties()->SetHeatScaleMidThreshold( (double)m_sliderHeatScaleMid->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}


void PanelVolume::OnSliderHeatScaleOffsetChanged( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      // double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      layer->GetProperties()->SetHeatScaleOffset( (double)m_sliderHeatScaleOffset->GetValue() / 100.0 * ( fMax + fMax ) - fMax );
    }
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
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->SetActiveFrame( event.GetInt() - 1 );
    }
  }
}

void PanelVolume::OnTextFrameChanged( wxCommandEvent& event )
{
  if ( m_listBoxLayers && m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    long value;
    if ( layer && m_textFrame->GetValue().ToLong( &value ) )
    {
      layer->SetActiveFrame( value - 1 );
    }
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
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->Lock( event.IsChecked() );
    }
  }
}

void PanelVolume::OnVolumeLockUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    event.Check( layer && layer->IsLocked() );
  }
}

void PanelVolume::OnVolumeCopySetting( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    m_layerCopied = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  }
}

void PanelVolume::OnVolumeCopySettingUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
}

void PanelVolume::OnVolumePasteSetting( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    layer->GetProperties()->CopySettings( m_layerCopied->GetProperties() );
  }
}

void PanelVolume::OnVolumePasteSettingUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    event.Enable( lc->Contains( m_layerCopied ) && m_layerCopied != layer );
  }
}


void PanelVolume::OnVolumePasteSettingAll( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    LayerMRI* layer = (LayerMRI*)lc->GetLayer( i );
    if ( layer->GetProperties()->GetColorMap() == m_layerCopied->GetProperties()->GetColorMap() )
      layer->GetProperties()->CopySettings( m_layerCopied->GetProperties() );
  }
}

void PanelVolume::OnVolumePasteSettingAllUpdateUI( wxUpdateUIEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  event.Enable( lc->Contains( m_layerCopied ) );
}

void PanelVolume::OnCheckContour( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    layer->GetProperties()->SetShowAsContour( event.IsChecked() );
  }
}


void PanelVolume::OnTextContourMin( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    // update both threshold input
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( !layer )
      return;
    
    double dMin = layer->GetProperties()->GetContourMinThreshold();
    double dMax = layer->GetProperties()->GetContourMaxThreshold();
    m_textContourMin->GetValue().ToDouble( &dMin );
    m_textContourMax->GetValue().ToDouble( &dMax );
    layer->GetProperties()->SetContourThreshold( dMin, dMax );
  }
}

void PanelVolume::OnTextContourMax( wxCommandEvent& event )
{
  /*
  double dvalue;
  if ( m_textContourMax->GetValue().ToDouble( &dvalue ) )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer && layer->GetProperties()->GetContourMaxThreshold() != dvalue )
    {
      layer->GetProperties()->SetContourMaxThreshold( dvalue );
    }
  }
  */
  OnTextContourMin( event );
}


void PanelVolume::OnSliderContourMin( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      layer->GetProperties()->SetContourMinThreshold( (double)m_sliderContourMin->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourMax( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      layer->GetProperties()->SetContourMaxThreshold( (double)m_sliderContourMax->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourMinChanging( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      UpdateTextValue( m_textContourMin, (double)m_sliderContourMin->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourMaxChanging( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      double fMin = layer->GetProperties()->GetMinValue();
      double fMax = layer->GetProperties()->GetMaxValue();
      UpdateTextValue( m_textContourMax, (double)m_sliderContourMax->GetValue() / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourSmooth( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetContourSmoothIterations( m_sliderContourSmooth->GetValue() );
    }
  }
}

void PanelVolume::OnSliderContourSmoothChanging( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      UpdateTextValue( m_textContourSmooth, m_sliderContourSmooth->GetValue() );
    }
  }  
}

void PanelVolume::OnTextContourSmooth( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    // update both threshold input
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( !layer )
      return;
    
    long n;
    if ( m_textContourSmooth->GetValue().ToLong( &n ) && n >= 0 && n <=50 )
      layer->GetProperties()->SetContourSmoothIterations( n );
  }
}

void PanelVolume::OnCheckDisplayVector( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetDisplayVector( event.IsChecked() );
      UpdateUI();
    } 
  }
}

void PanelVolume::OnCheckDisplayTensor( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetDisplayTensor( event.IsChecked() );
      UpdateUI();
    } 
  }
}

void PanelVolume::OnChoiceInversion( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
}

void PanelVolume::OnChoiceRepresentation( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
}

void PanelVolume::OnCheckHeatScaleClearHigh( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleClearHigh( event.IsChecked() );
    }
  }
}

void PanelVolume::OnCheckHeatScaleTruncate( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleTruncate( event.IsChecked() );
    }
  }
}

void PanelVolume::OnCheckHeatScaleInvert( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetHeatScaleInvert( event.IsChecked() );
    }
  }
}

void PanelVolume::OnCheckUseImageColorMap( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetContourUseImageColorMap( event.IsChecked() );
    }
  }
}

void PanelVolume::OnCheckContourExtractAll( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetContourExtractAllRegions( event.IsChecked() );
    }
  }
}

void PanelVolume::OnColorContour( wxColourPickerEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      wxColour c = event.GetColour();
      layer->GetProperties()->SetContourColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
    }
  }
}

/*
void PanelVolume::OnColorSurface( wxColourPickerEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
  if ( layer )
  {
    wxColour c = event.GetColour();
    layer->GetSurfaceProperties()->SetColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}
*/
void PanelVolume::OnCheckHideInfo( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetShowInfo( !event.IsChecked() );
    }
  }
}

void PanelVolume::OnCheckShowLabelOutline( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetShowLabelOutline( event.IsChecked() );
    }
  }
}

void PanelVolume::OnChoiceUpSampleMethod( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      layer->GetProperties()->SetUpSampleMethod( event.GetSelection() );
    }
  }
}

void PanelVolume::OnButtonSaveContour( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      wxFileDialog dlg( this, _("Save iso-surface"), 
                        MainWindow::GetMainWindowPointer()->AutoSelectLastDir( _("mri") ),
                        wxString( layer->GetName() ) + _(".vtk"),
                        _("VTK files (*.vtk)|*.vtk|All files (*.*)|*.*"),
                        wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
      if ( dlg.ShowModal() == wxID_OK )
      {
        if ( !layer->SaveContourToFile( dlg.GetPath().char_str() ) )
        {
          wxMessageDialog mdlg( this, _("Can not save surface to file."), 
                               _("Error"), wxOK );
          mdlg.ShowModal();
        }
      }
    }
  }  
}

