/**
 * @file  PanelSceneSetting.h
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
#include "PanelSceneSetting.h"
#include <wx/xrc/xmlres.h>
#include "MainWindow.h"
#include "LayerCollection.h"
#include "Layer.h"
#include "LayerROI.h"
#include "LayerPropertiesROI.h"

BEGIN_EVENT_TABLE( PanelSceneSetting, wxPanel )
END_EVENT_TABLE()


PanelSceneSetting::PanelSceneSetting( wxWindow* parent ) :
    Listener( "PanelSceneSetting" ),
    Broadcaster( "PanelSceneSetting" ),
    m_bUINeedUpdate( false )
{
  wxXmlResource::Get()->LoadPanel( this, parent, wxT("ID_PANEL_SCENE_SETTINGS") );

  UpdateUI();
}

PanelSceneSetting::~PanelSceneSetting()
{}

void PanelSceneSetting::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// MainWindow* mainwnd = MainWindow::GetMainWindow();
// LayerCollection* lc = mainwnd->GetLayerCollection();
  /* if ( iMsg == "LayerAdded" )
   {
    Layer* layer = ( Layer* )iData;
    if ( layer && layer->IsTypeOf( "ROI" ) )
    {
     m_listBoxLayers->Insert( layer->GetName(), 0, (void*)layer );
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
   */
}

void PanelSceneSetting::OnSliderX( wxScrollEvent& event )
{}

void PanelSceneSetting::OnSliderY( wxScrollEvent& event )
{}

void PanelSceneSetting::OnSliderZ( wxScrollEvent& event )
{}

void PanelSceneSetting::UpdateUI( bool bForce )
{
  if ( bForce )
    DoUpdateUI();
  else
    m_bUINeedUpdate = true;
}

void PanelSceneSetting::DoUpdateUI()
{
  /*
  bool bHasROI = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
  bool bHasVolume = !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty();
  LayerROI* layer = NULL;
  if ( bHasROI )
  {
   layer = ( LayerROI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
   if ( layer )
   {
    m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
    double* rgb = layer->GetProperties()->GetColor();
    m_colorPicker->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    m_textFileName->ChangeValue( layer->GetFileName() );
    m_textFileName->SetInsertionPointEnd();
    m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
   }
   
   LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
   lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
  }
  MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  m_btnNew->Enable( bHasVolume );
  m_btnLoad->Enable( bHasVolume );
  m_btnDelete->Enable( bHasROI && !mainWnd->IsSaving() ); 
  m_btnMoveUp->Enable( bHasROI && m_listBoxLayers->GetSelection() != 0 );
  m_btnMoveDown->Enable( bHasROI && m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
  m_btnSave->Enable( bHasROI && layer && layer->IsModified() && !mainWnd->IsSaving() );
  */
}


void PanelSceneSetting::OnInternalIdle()
{
  if ( m_bUINeedUpdate )
  {
    DoUpdateUI();
    m_bUINeedUpdate = false;
  }
  wxPanel::OnInternalIdle();
}
