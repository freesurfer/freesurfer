/**
 * @file  PanelVolume.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/05/20 16:28:32 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
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
/*	EVT_IDLE(PanelVolume::OnIdle)
	EVT_CHECKBOX		(XRCID(wxT("ID_CHECKBOX_SLICE_Y")),			PanelVolume::OnCheckBoxSlice)
	EVT_CHECKBOX		(XRCID(wxT("ID_CHECKBOX_SLICE_Z")),			PanelVolume::OnCheckBoxSlice)
	EVT_TEXT			(XRCID(wxT("ID_TEXTCTRL_SLICE_X")),			PanelVolume::OnEditSlice) 
	EVT_TEXT			(XRCID(wxT("ID_TEXTCTRL_SLICE_Y")),			PanelVolume::OnEditSlice) 
	EVT_TEXT			(XRCID(wxT("ID_TEXTCTRL_SLICE_Z")),			PanelVolume::OnEditSlice) 
	EVT_COMMAND_SCROLL	(XRCID(wxT("ID_SLIDER_SLICE_X")),			PanelVolume::OnSliderSlice)
	EVT_COMMAND_SCROLL	(XRCID(wxT("ID_SLIDER_SLICE_Y")),			PanelVolume::OnSliderSlice)*/
	EVT_LISTBOX			( XRCID( wxT( "ID_LISTBOX_VOLUMES" ) ),			PanelVolume::OnLayerSelectionChanged )	
	EVT_CHECKLISTBOX	( XRCID( wxT( "ID_LISTBOX_VOLUMES" ) ),			PanelVolume::OnLayerVisibilityChanged )
	EVT_LISTBOX_DCLICK 	( XRCID( wxT( "ID_LISTBOX_VOLUMES" ) ),			PanelVolume::OnListDoubleClicked )
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_OPACITY" ) ),			PanelVolume::OnSliderOpacityChanged )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_LOAD" ) ), 			PanelVolume::OnButtonLoad )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_MOVE_UP" ) ),			PanelVolume::OnButtonMoveUp )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_MOVE_DOWN" ) ),		PanelVolume::OnButtonMoveDown )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_DELETE" ) ),			PanelVolume::OnButtonDelete )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_NEW" ) ),				PanelVolume::OnButtonNew )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_LOAD" ) ),				PanelVolume::OnButtonLoad )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SAVE" ) ),				PanelVolume::OnButtonSave )
	EVT_CHECKBOX		( XRCID( wxT( "ID_CHECKBOX_CLEAR_BACKGROUND" ) ),	PanelVolume::OnCheckClearBackground )
	EVT_CHECKBOX		( XRCID( wxT( "ID_CHECKBOX_SMOOTH" ) ),			PanelVolume::OnCheckSmooth )
	EVT_CHOICE			( XRCID( wxT( "ID_CHOICE_COLORMAP" ) ),			PanelVolume::OnChoiceColorMap )
	EVT_CHOICE			( XRCID( wxT( "ID_CHOICE_LUT" ) ),				PanelVolume::OnChoiceLUT )	
	EVT_CHOICE			( XRCID( wxT( "ID_CHOICE_DIRECTION_CODE" ) ),	PanelVolume::OnChoiceDirectionCode )	
	EVT_LISTBOX			( XRCID( wxT( "ID_LISTBOX_COLORTABLE" ) ),		PanelVolume::OnColorSelectionChanged )
	EVT_LISTBOX_DCLICK	( XRCID( wxT( "ID_LISTBOX_COLORTABLE" ) ),		PanelVolume::OnColorSelectionChanged )	
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_DRAW_VALUE" ) ),			PanelVolume::OnTextFillValueChanged )	
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_WINDOW" ) ),				PanelVolume::OnTextWindowChanged )	
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_LEVEL" ) ),				PanelVolume::OnTextLevelChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_HEATSCALE_MIN" ) ),		PanelVolume::OnTextHeatScaleMinChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_HEATSCALE_MID" ) ),		PanelVolume::OnTextHeatScaleMidChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_HEATSCALE_MAX" ) ),		PanelVolume::OnTextHeatScaleMaxChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_JETSCALE_MIN" ) ),		PanelVolume::OnTextMinJetScaleChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_JETSCALE_MAX" ) ),		PanelVolume::OnTextMaxJetScaleChanged )
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_OPACITY" ) ),			PanelVolume::OnTextOpacityChanged )
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_WINDOW" ) ),			PanelVolume::OnSliderWindowChanged )
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_LEVEL" ) ),			PanelVolume::OnSliderLevelChanged )	
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_HEATSCALE_MIN" ) ),	PanelVolume::OnSliderHeatScaleMinChanged )	
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_HEATSCALE_MID" ) ),	PanelVolume::OnSliderHeatScaleMidChanged )	
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_HEATSCALE_MAX" ) ),	PanelVolume::OnSliderHeatScaleMaxChanged )	
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_JETSCALE_MIN" ) ),		PanelVolume::OnSliderMinJetScaleChanged )	
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_JETSCALE_MAX" ) ),		PanelVolume::OnSliderMaxJetScaleChanged )		
	EVT_TEXT			( XRCID( wxT( "ID_TEXT_FRAME" ) ),				PanelVolume::OnTextFrameChanged )
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_FRAME" ) ),			PanelVolume::OnSliderFrameChanged )	
END_EVENT_TABLE()


PanelVolume::PanelVolume( wxWindow* parent ) : Listener( "PanelVolume" ), Broadcaster( "PanelVolume" )
{
	m_curCTAB = NULL;
	m_bUINeedUpdate = false;
	
	wxXmlResource::Get()->LoadPanel( this, parent, wxT("ID_PANEL_VOLUME") );
	m_btnNew = XRCCTRL( *this, "ID_BUTTON_NEW", wxButton );
	m_btnDelete = XRCCTRL( *this, "ID_BUTTON_DELETE", wxButton );
	m_btnMoveUp = XRCCTRL( *this, "ID_BUTTON_MOVE_UP", wxButton );
	m_btnMoveDown = XRCCTRL( *this, "ID_BUTTON_MOVE_DOWN", wxButton );
	m_btnSave = XRCCTRL( *this, "ID_BUTTON_SAVE", wxButton );
	m_listBoxLayers = XRCCTRL( *this, "ID_LISTBOX_VOLUMES", wxCheckListBox );
	m_sliderOpacity = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
	m_textOpacity = XRCCTRL( *this, "ID_TEXT_OPACITY", wxTextCtrl );
	m_checkClearBackground = XRCCTRL( *this, "ID_CHECKBOX_CLEAR_BACKGROUND", wxCheckBox );
	m_checkSmooth = XRCCTRL( *this, "ID_CHECKBOX_SMOOTH", wxCheckBox );
	m_listColorTable = XRCCTRL( *this, "ID_LISTBOX_COLORTABLE", wxListBox );
	m_choiceColorMap = XRCCTRL( *this, "ID_CHOICE_COLORMAP", wxChoice );
	m_choiceLUT = XRCCTRL( *this, "ID_CHOICE_LUT", wxChoice );
	m_textDrawValue = XRCCTRL( *this, "ID_TEXT_DRAW_VALUE", wxTextCtrl );
	m_sliderWindow = XRCCTRL( *this, "ID_SLIDER_WINDOW", wxSlider );
	m_sliderLevel = XRCCTRL( *this, "ID_SLIDER_LEVEL", wxSlider );
	m_textWindow = XRCCTRL( *this, "ID_TEXT_WINDOW", wxTextCtrl );
	m_textLevel = XRCCTRL( *this, "ID_TEXT_LEVEL", wxTextCtrl );
	m_textFileName = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
	m_textHeatScaleMin = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MIN", wxTextCtrl );
	m_textHeatScaleMid = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MID", wxTextCtrl );
	m_textHeatScaleMax = XRCCTRL( *this, "ID_TEXT_HEATSCALE_MAX", wxTextCtrl );
	m_sliderHeatScaleMin = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MIN", wxSlider );
	m_sliderHeatScaleMid = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MID", wxSlider );
	m_sliderHeatScaleMax = XRCCTRL( *this, "ID_SLIDER_HEATSCALE_MAX", wxSlider );
	m_textJetScaleMin = XRCCTRL( *this, "ID_TEXT_JETSCALE_MIN", wxTextCtrl );
	m_textJetScaleMax = XRCCTRL( *this, "ID_TEXT_JETSCALE_MAX", wxTextCtrl );
	m_sliderJetScaleMin = XRCCTRL( *this, "ID_SLIDER_JETSCALE_MIN", wxSlider );
	m_sliderJetScaleMax = XRCCTRL( *this, "ID_SLIDER_JETSCALE_MAX", wxSlider );
	m_choiceDirectionCode = XRCCTRL( *this, "ID_CHOICE_DIRECTION_CODE", wxChoice );
	m_colorIndicator = XRCCTRL( *this, "ID_COLOR_INDICATOR", wxColorIndicator );
	m_textFrame = XRCCTRL( *this, "ID_TEXT_FRAME", wxTextCtrl );
	m_sliderFrame = XRCCTRL( *this, "ID_SLIDER_FRAME", wxSlider );
	
	m_luts = MainWindow::GetMainWindowPointer()->GetLUTData();
	for ( int i = 0; i < m_luts->GetCount(); i++ )
	{
		m_choiceLUT->Append( m_luts->GetName( i ) );
	}
//	m_choiceLUT->Append( "Edit/add lookup table..." );
	
	m_choiceDirectionCode->Append( "RAS -> RGB" );
	m_choiceDirectionCode->Append( "RAS -> RBG" );
	m_choiceDirectionCode->Append( "RAS -> GRB" );
	m_choiceDirectionCode->Append( "RAS -> GBR" );
	m_choiceDirectionCode->Append( "RAS -> BRG" );
	m_choiceDirectionCode->Append( "RAS -> BGR" );
	
	m_widgetlistGrayScale.push_back( m_checkClearBackground );
	m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_WINDOW", wxStaticText ) );
	m_widgetlistGrayScale.push_back( XRCCTRL( *this, "ID_STATIC_LEVEL", wxStaticText ) );
	m_widgetlistGrayScale.push_back( m_textWindow );
	m_widgetlistGrayScale.push_back( m_textLevel );
	m_widgetlistGrayScale.push_back( m_sliderWindow );
	m_widgetlistGrayScale.push_back( m_sliderLevel );
	
	m_widgetlistHeatScale.push_back( m_textHeatScaleMin );
	m_widgetlistHeatScale.push_back( m_textHeatScaleMid );
	m_widgetlistHeatScale.push_back( m_textHeatScaleMax );
	m_widgetlistHeatScale.push_back( m_sliderHeatScaleMin );
	m_widgetlistHeatScale.push_back( m_sliderHeatScaleMid );
	m_widgetlistHeatScale.push_back( m_sliderHeatScaleMax );
	m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MIN", wxStaticText ) );
	m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MID", wxStaticText ) );
	m_widgetlistHeatScale.push_back( XRCCTRL( *this, "ID_STATIC_HEATSCALE_MAX", wxStaticText ) );
	
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
	
	MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->AddListener( this );
	UpdateUI();
}

PanelVolume::~PanelVolume()
{
}

void PanelVolume::DoListenToMessage( std::string const iMsg, void* iData )
{
//	MainWindow* mainwnd = MainWindow::GetMainWindow();
//	LayerCollection* lc = mainwnd->GetLayerCollection();
	if ( iMsg == "LayerAdded" )
	{
		Layer* layer = ( Layer* )iData;
		if ( layer && layer->IsTypeOf( "MRI" ) )
		{
			m_listBoxLayers->Insert( layer->GetName(), 0, (void*)layer );
			m_listBoxLayers->Check( 0 );
			m_listBoxLayers->SetSelection( 0 );
			UpdateUI();
		}
	}
	else if ( iMsg == "LayerModified" || iMsg == /*"WindowLevelChanged"*/ "LayerActorUpdated" )
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
		layer->SetVisible( m_listBoxLayers->IsChecked( nItem ) );
	}
}

void PanelVolume::OnButtonMoveUp( wxCommandEvent& event )
{
	LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
	int nSel = m_listBoxLayers->GetSelection();
	if ( lc && nSel != wxNOT_FOUND )
	{	
		Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );
		
		if ( nSel != 0 )
		{
			void* old_ptr = m_listBoxLayers->GetClientData( nSel - 1 );
			wxString old_name = m_listBoxLayers->GetString( nSel - 1 );
			bool old_checked = m_listBoxLayers->IsChecked( nSel - 1 );
			
			m_listBoxLayers->Delete( nSel - 1 );
			m_listBoxLayers->Insert( old_name, nSel, old_ptr );
			m_listBoxLayers->Check( nSel, old_checked ); 
			
			if ( layer )
				lc->MoveLayerUp( layer );
			
			UpdateUI();
		}
	}
}

void PanelVolume::OnButtonMoveDown( wxCommandEvent& event )
{
	LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
	int nSel = m_listBoxLayers->GetSelection();
	if ( lc && nSel != wxNOT_FOUND )
	{
		Layer* layer = ( Layer* )( void* )m_listBoxLayers->GetClientData( nSel );
		
		if ( nSel != (int)m_listBoxLayers->GetCount()-1 )
		{
			void* old_ptr = m_listBoxLayers->GetClientData( nSel + 1 );
			wxString old_name = m_listBoxLayers->GetString( nSel + 1 );
			bool old_checked = m_listBoxLayers->IsChecked( nSel + 1 );
			
			m_listBoxLayers->Delete( nSel + 1 );
			m_listBoxLayers->Insert( old_name, nSel, old_ptr ); 
			m_listBoxLayers->Check( nSel, old_checked ); 
			
			if ( layer )
				lc->MoveLayerDown( layer );
			
			UpdateUI();
		}
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
			wxString msg = "Volume has been modified. Do you want to close it without saving?";
			wxMessageDialog dlg( this, msg, "Close", wxYES_NO | wxICON_QUESTION | wxNO_DEFAULT );
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
			wxTextEntryDialog dlg( this, "Enter the name of the volume:", "Volume", layer->GetName() );
			if ( dlg.ShowModal() == wxID_OK )
			{
				layer->SetName( dlg.GetValue().Trim( true ).Trim( false ).c_str() );
				bool bChecked = m_listBoxLayers->IsChecked( event.GetInt() );
				m_listBoxLayers->SetString( event.GetInt(), layer->GetName() );
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
		// to add/edit lookup table
		DialogEditLookupTable dlg( this );
		dlg.ShowModal();
		return;
	}	
	
	if ( m_listBoxLayers->GetSelection() != wxNOT_FOUND )
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
				m_listColorTable->Append( wxString::Format( "%d: %s", i, name ) );
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
	for ( int i = 0; i < (int)list.size(); i++ )
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
	LayerMRI* layer = NULL;
	LayerPropertiesMRI::ColorMapType nColorMap = LayerPropertiesMRI::NoColorMap;
	if ( bHasVolume )
	{
		layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
		if ( layer )
		{
			m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
			UpdateTextValue( m_textOpacity, layer->GetProperties()->GetOpacity() );
			m_checkClearBackground->SetValue( layer->GetProperties()->GetClearZero() );
			m_checkSmooth->SetValue( layer->GetProperties()->GetTextureSmoothing() );
			
			m_choiceColorMap->SetSelection( layer->GetProperties()->GetColorMap() );
			m_choiceLUT->SetSelection( m_luts->GetIndex( layer->GetProperties()->GetLUTCTAB() ) );
	
			UpdateTextValue( m_textDrawValue, (double)layer->GetFillValue() );	
			UpdateTextValue( m_textWindow, layer->GetProperties()->GetWindow() );
			UpdateTextValue( m_textLevel, layer->GetProperties()->GetLevel() );
			double* windowrange = layer->GetProperties()->GetWindowRange();
			double* levelrange = layer->GetProperties()->GetLevelRange();
			m_sliderWindow->SetValue( (int)( ( layer->GetProperties()->GetWindow() - windowrange[0] ) / ( windowrange[1] - windowrange[0] ) * 100 ) );
			m_sliderLevel->SetValue( (int)( ( layer->GetProperties()->GetLevel() - levelrange[0] ) / ( levelrange[1] - levelrange[0] ) * 100 ) );
			
			if ( layer->IsTypeOf( "DTI" ) )
				m_textFileName->ChangeValue( ((LayerDTI*)layer)->GetVectorFileName() );
			else
				m_textFileName->ChangeValue( layer->GetFileName() );
			m_textFileName->SetInsertionPointEnd();
			m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
			
			double fMin = layer->GetProperties()->GetMinValue();
			double fMax = layer->GetProperties()->GetMaxValue();
			m_sliderHeatScaleMin->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMinThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderHeatScaleMid->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMidThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			m_sliderHeatScaleMax->SetValue( (int)( ( layer->GetProperties()->GetHeatScaleMaxThreshold() - fMin ) / ( fMax - fMin ) * 100 ) );
			UpdateTextValue( m_textHeatScaleMin, layer->GetProperties()->GetHeatScaleMinThreshold() );
			UpdateTextValue( m_textHeatScaleMid, layer->GetProperties()->GetHeatScaleMidThreshold() );
			UpdateTextValue( m_textHeatScaleMax, layer->GetProperties()->GetHeatScaleMaxThreshold() );
			
			double fJetMin = fMin - (fMax-fMin)/4;
			double fJetMax = fMax + (fMax-fMin)/4;
			m_sliderJetScaleMin->SetValue( (int)( ( layer->GetProperties()->GetMinJetScaleWindow() - fJetMin ) / ( fJetMax - fJetMin ) * 100 ) );
			m_sliderJetScaleMax->SetValue( (int)( ( layer->GetProperties()->GetMaxJetScaleWindow() - fJetMin ) / ( fJetMax - fJetMin ) * 100 ) );
			UpdateTextValue( m_textJetScaleMin, layer->GetProperties()->GetMinJetScaleWindow() );
			UpdateTextValue( m_textJetScaleMax, layer->GetProperties()->GetMaxJetScaleWindow() );
					
			m_choiceColorMap->Clear();
			m_choiceColorMap->Append( "Grayscale", (void*)LayerPropertiesMRI::Grayscale );
			m_choiceColorMap->Append( "Heat", (void*)LayerPropertiesMRI::Heat );
			m_choiceColorMap->Append( "Jet", (void*)LayerPropertiesMRI::Jet );
			if ( layer->IsTypeOf( "DTI" ) )
			{
				m_choiceColorMap->Append( "Direction-coded", (void*)LayerPropertiesMRI::DirectionCoded );
				m_choiceDirectionCode->SetSelection( ((LayerDTI*)layer)->GetProperties()->GetDirectionCode() );
			}
			else
			{
				m_choiceColorMap->Append( "Lookup Table", (void*)LayerPropertiesMRI::LUT );		
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
		}
	}
	MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
	m_btnNew->Enable( bHasVolume );
	m_btnDelete->Enable( bHasVolume && !mainWnd->IsSaving() );	
	m_btnMoveUp->Enable( bHasVolume && m_listBoxLayers->GetSelection() != 0 );
	m_btnMoveDown->Enable( bHasVolume && m_listBoxLayers->GetSelection() != ( (int)m_listBoxLayers->GetCount() - 1 ) );
	m_btnSave->Enable( bHasVolume && layer && layer->IsModified() && !mainWnd->IsSaving() );	
	
	ShowWidgets( m_widgetlistGrayScale, bHasVolume && nColorMap == LayerPropertiesMRI::Grayscale );
	ShowWidgets( m_widgetlistHeatScale, bHasVolume && nColorMap == LayerPropertiesMRI::Heat );
	ShowWidgets( m_widgetlistJetScale, bHasVolume && nColorMap == LayerPropertiesMRI::Jet );
	ShowWidgets( m_widgetlistLUT, bHasVolume && nColorMap == LayerPropertiesMRI::LUT );
	ShowWidgets( m_widgetlistDirectionCode, bHasVolume && nColorMap == LayerPropertiesMRI::DirectionCoded );
	ShowWidgets( m_widgetlistFrame, bHasVolume && layer && layer->GetNumberOfFrames() > 1 );
	
	Layout();
	if ( layer && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
	{
		PopulateColorTable( layer->GetProperties()->GetLUTCTAB() );
		
		UpdateColorIndicator();
	}
}

void PanelVolume::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
	wxString value_strg = ( (wxString)"" << dvalue );
	if ( value_strg != ctrl->GetValue() && (value_strg + ".") != ctrl->GetValue() )
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
	if ( strg.Left( strg.Find( ":" ) ).ToDouble( &dValue ) )
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
			if ( strg.Left( strg.Find( ":" ) ).ToDouble( &dvalue ) && dvalue == (double)nvalue )
			{
				m_listColorTable->SetFirstItem( i );
				m_listColorTable->SetSelection( i );
				break;
			}
		}
		
		UpdateColorIndicator();
				
		LayerMRI* layer = ( LayerMRI* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
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
		//	this->BlockListen( true );
			layer->GetProperties()->SetWindow( dvalue );
		//	this->BlockListen( false );		
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


