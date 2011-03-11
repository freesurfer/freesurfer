/**
 * @file  PanelROI.h
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
#include "PanelROI.h"
#include <wx/xrc/xmlres.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerROI.h"
#include "LayerPropertiesROI.h"

BEGIN_EVENT_TABLE( PanelROI, wxPanel )
  EVT_LISTBOX         ( XRCID( "ID_LISTBOX_ROIS" ),     PanelROI::OnLayerSelectionChanged )
  EVT_CHECKLISTBOX    ( XRCID( "ID_LISTBOX_ROIS" ),     PanelROI::OnLayerVisibilityChanged )
  EVT_COMMAND_SCROLL  ( XRCID( "ID_SLIDER_OPACITY" ),   PanelROI::OnSliderOpacity )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_MOVE_UP" ),   PanelROI::OnButtonMoveUp )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_MOVE_DOWN" ), PanelROI::OnButtonMoveDown )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_DELETE" ),    PanelROI::OnButtonDelete )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_NEW" ),       PanelROI::OnButtonNew )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_LOAD" ),      PanelROI::OnButtonLoad )
  EVT_BUTTON          ( XRCID( "ID_BUTTON_SAVE" ),      PanelROI::OnButtonSave )
  EVT_MENU            ( XRCID( "ID_ROI_CLOSE" ),        PanelROI::OnButtonDelete )
  EVT_UPDATE_UI       ( XRCID( "ID_ROI_CLOSE" ),        PanelROI::OnROICloseUpdateUI )
  EVT_MENU            ( XRCID( "ID_ROI_MOVE_UP" ),      PanelROI::OnButtonMoveUp )
  EVT_UPDATE_UI       ( XRCID( "ID_ROI_MOVE_UP" ),      PanelROI::OnMoveUpUpdateUI )
  EVT_MENU            ( XRCID( "ID_ROI_MOVE_DOWN" ),    PanelROI::OnButtonMoveDown )
  EVT_UPDATE_UI       ( XRCID( "ID_ROI_MOVE_DOWN" ),    PanelROI::OnMoveDownUpdateUI )
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLOR_PICKER" ),  PanelROI::OnColorChanged )
END_EVENT_TABLE()


PanelROI::PanelROI( wxWindow* parent ) :
    Listener( "PanelROI" ),
    Broadcaster( "PanelROI" ),
    m_bUINeedUpdate( false )
{
  wxXmlResource::Get()->LoadPanel( this, parent, wxT("ID_PANEL_ROI") );
  m_listBoxLayers = XRCCTRL( *this, "ID_LISTBOX_ROIS", wxCheckListBox );
  m_sliderOpacity = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
  m_btnNew = XRCCTRL( *this, "ID_BUTTON_NEW", wxButton );
  m_btnLoad = XRCCTRL( *this, "ID_BUTTON_LOAD", wxButton );
  m_btnSave = XRCCTRL( *this, "ID_BUTTON_SAVE", wxButton );
  m_btnDelete = XRCCTRL( *this, "ID_BUTTON_DELETE", wxButton );
  m_btnMoveUp = XRCCTRL( *this, "ID_BUTTON_MOVE_UP", wxButton );
  m_btnMoveDown = XRCCTRL( *this, "ID_BUTTON_MOVE_DOWN", wxButton );
  m_colorPicker = XRCCTRL( *this, "ID_COLOR_PICKER", wxColourPickerCtrl );
  m_textFileName = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" )->AddListener( this );
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->AddListener( this );

  UpdateUI();
}

PanelROI::~PanelROI()
{}

void PanelROI::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// MainWindow* mainwnd = MainWindow::GetMainWindow();
// LayerCollection* lc = mainwnd->GetLayerCollection();
  if ( iMsg == "LayerAdded" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "ROI" ) )
    {
      m_listBoxLayers->Insert( wxString::FromAscii( layer->GetName() ), 0, (void*)layer );
      m_listBoxLayers->Check( 0 );
      m_listBoxLayers->SetSelection( 0 );
    }

    UpdateUI();
  }
  else if ( iMsg == "LayerMoved" )
  {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "ROI" ) )
    {
      UpdateLayerList( layer );
    }
  }
  else if ( iMsg == "LayerModified" )
  {
    UpdateUI();
  }
}

void PanelROI::OnSliderOpacity( wxScrollEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerROI* layer = ( LayerROI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
      layer->GetProperties()->SetOpacity( event.GetPosition() / 100.0 );
  }
}

void PanelROI::OnLayerSelectionChanged( wxCommandEvent& event )
{
  if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  UpdateUI();
}

void PanelROI::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void PanelROI::DoUpdateUI()
{
  bool bHasROI = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  bool bHasVolume = !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty();
  wxWindowList children = GetChildren();
  wxWindowList::iterator it = children.begin(), end = children.end();
  for (; it != end; it++)
  {
    if ( !(*it)->IsKindOf(CLASSINFO(wxToolBar) ) && *it != m_listBoxLayers )
      (*it)->Enable( bHasROI );
  }

  LayerROI* layer = NULL;
  if ( bHasROI )
  {
    layer = ( LayerROI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
    if ( layer )
    {
      m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
      double* rgb = layer->GetProperties()->GetColor();
      m_colorPicker->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
      m_textFileName->ChangeValue( wxString::FromAscii( layer->GetFileName() ) );
      m_textFileName->SetInsertionPointEnd();
      m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
    }

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
    lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  m_btnNew->Enable( bHasVolume );
  m_btnLoad->Enable( bHasVolume );
  m_btnDelete->Enable( bHasROI && !mainWnd->IsProcessing() );
  m_btnMoveUp->Enable( bHasROI && m_listBoxLayers->GetSelection() != 0 );
  m_btnMoveDown->Enable( bHasROI && m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
  m_btnSave->Enable( bHasROI && layer && layer->IsModified() && !mainWnd->IsProcessing() );
}

void PanelROI::OnLayerVisibilityChanged( wxCommandEvent& event )
{
  int nItem = event.GetInt();
  Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nItem );
  if ( layer )
  {
    layer->SetVisible( m_listBoxLayers->IsChecked( nItem ) );
  }
}

void PanelROI::OnButtonMoveUp( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

    if ( layer )
      lc->MoveLayerUp( layer );
  }
}

void PanelROI::OnButtonMoveDown( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );

    if ( layer )
      lc->MoveLayerDown( layer );
  }
}


void PanelROI::UpdateLayerList( Layer* layer )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
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

void PanelROI::OnButtonDelete( wxCommandEvent& event )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
  int nSel = m_listBoxLayers->GetSelection();
  if ( lc && nSel != wxNOT_FOUND )
  {
    Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );
    if ( ((LayerROI*)layer)->IsModified() )
    {
      wxString msg = _("ROI has been modified. Do you want to close it without saving?");
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

void PanelROI::OnButtonNew( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->NewROI();
}

void PanelROI::OnButtonLoad( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->LoadROI();
}

void PanelROI::OnButtonSave( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveROI();
}

void PanelROI::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxPanel::OnInternalIdle();
}

void PanelROI::OnColorChanged( wxColourPickerEvent& event )
{
  LayerROI* roi = ( LayerROI* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" )->GetActiveLayer();
  if ( roi )
  {
    wxColour c = event.GetColour();
    roi->GetProperties()->SetColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

void PanelROI::OnROICloseUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND &&
                !MainWindow::GetMainWindowPointer()->IsProcessing() );
}

void PanelROI::OnMoveUpUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND &&
                m_listBoxLayers->GetSelection() != 0 );
}

void PanelROI::OnMoveDownUpdateUI( wxUpdateUIEvent& event )
{
  event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND &&
                m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
}
