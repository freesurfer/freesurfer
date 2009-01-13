/**
 * @file  PanelSurface.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/01/13 21:19:34 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2009,
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

BEGIN_EVENT_TABLE( PanelSurface, wxPanel )
	EVT_LISTBOX			( XRCID( wxT( "ID_LISTBOX_SURFACE" ) ),			PanelSurface::OnLayerSelectionChanged )	
	EVT_CHECKLISTBOX	( XRCID( wxT( "ID_LISTBOX_SURFACE" ) ),			PanelSurface::OnLayerVisibilityChanged )
	EVT_COMMAND_SCROLL	( XRCID( wxT( "ID_SLIDER_OPACITY" ) ),			PanelSurface::OnSliderOpacity )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_LOAD" ) ),				PanelSurface::OnButtonLoad )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_DELETE" ) ),			PanelSurface::OnButtonDelete )
	EVT_MENU			( XRCID( wxT( "ID_SURFACE_CLOSE" ) ),			PanelSurface::OnButtonDelete )
	EVT_UPDATE_UI		( XRCID( wxT( "ID_SURFACE_CLOSE" ) ),			PanelSurface::OnSurfaceCloseUpdateUI )
	
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SURFACE_MAIN" ) ),		PanelSurface::OnButtonSurfaceMain )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SURFACE_INFLATED" ) ),	PanelSurface::OnButtonSurfaceInflated )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SURFACE_WHITE" ) ),	PanelSurface::OnButtonSurfaceWhite )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SURFACE_PIAL" ) ),		PanelSurface::OnButtonSurfacePial )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_SURFACE_ORIGINAL" ) ),	PanelSurface::OnButtonSurfaceOriginal )
	EVT_COLOURPICKER_CHANGED	( XRCID( wxT( "ID_COLOR_PICKER" ) ), 	PanelSurface::OnColorChanged )
	EVT_COLOURPICKER_CHANGED	( XRCID( wxT( "ID_COLOR_PICKER_EDGE" ) ), 	PanelSurface::OnEdgeColorChanged )
	EVT_SPINCTRL		( XRCID( wxT( "ID_SPIN_EDGE_THICKNESS" ) ),		PanelSurface::OnSpinEdgeThickness )
	EVT_CHOICE			( XRCID( wxT( "ID_CHOICE_VECTORS" ) ), 			PanelSurface::OnChoiceVector )
END_EVENT_TABLE()


PanelSurface::PanelSurface( wxWindow* parent ) : 
		Listener( "PanelSurface" ), 
		Broadcaster( "PanelSurface" ),
		m_bUINeedUpdate( false )
{
	wxXmlResource::Get()->LoadPanel( this, parent, wxT("ID_PANEL_SURFACE") );
	m_listBoxLayers = XRCCTRL( *this, "ID_LISTBOX_SURFACE", wxCheckListBox );
	m_sliderOpacity = XRCCTRL( *this, "ID_SLIDER_OPACITY", wxSlider );
	m_btnNew 		= XRCCTRL( *this, "ID_BUTTON_NEW", wxButton );
	m_btnLoad 		= XRCCTRL( *this, "ID_BUTTON_LOAD", wxButton );
	m_btnSave 		= XRCCTRL( *this, "ID_BUTTON_SAVE", wxButton );
	m_btnDelete 	= XRCCTRL( *this, "ID_BUTTON_DELETE", wxButton );
	
	m_btnSurfaceMain 		= XRCCTRL( *this, "ID_BUTTON_SURFACE_MAIN", wxButton );
	m_btnSurfaceInflated 	= XRCCTRL( *this, "ID_BUTTON_SURFACE_INFLATED", wxButton );
	m_btnSurfaceWhite 		= XRCCTRL( *this, "ID_BUTTON_SURFACE_WHITE", wxButton );
	m_btnSurfacePial 		= XRCCTRL( *this, "ID_BUTTON_SURFACE_PIAL", wxButton );
	m_btnSurfaceOriginal 	= XRCCTRL( *this, "ID_BUTTON_SURFACE_ORIGINAL", wxButton );
	
	m_colorPicker 			= XRCCTRL( *this, "ID_COLOR_PICKER", wxColourPickerCtrl );
	m_colorPickerEdge 		= XRCCTRL( *this, "ID_COLOR_PICKER_EDGE", wxColourPickerCtrl );
	m_textFileName 			= XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
	m_spinEdgeThickness		= XRCCTRL( *this, "ID_SPIN_EDGE_THICKNESS", wxSpinCtrl );
	
	m_choiceVector			= XRCCTRL( *this, "ID_CHOICE_VECTORS", wxChoice );
	
	MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->AddListener( this );
	
	UpdateUI();
}

PanelSurface::~PanelSurface()
{
}

void PanelSurface::DoListenToMessage( std::string const iMsg, void* iData )
{
//	MainWindow* mainwnd = MainWindow::GetMainWindow();
//	LayerCollection* lc = mainwnd->GetLayerCollection();
	if ( iMsg == "LayerAdded" )
	{
		Layer* layer = ( Layer* )iData;
		if ( layer && layer->IsTypeOf( "Surface" ) )
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
		if ( layer && layer->IsTypeOf( "Surface" ) )
		{
			UpdateLayerList( layer );
		}
	}
	else if ( iMsg == "LayerModified" )
	{
		UpdateUI();
	}
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
	bool bHasSurface = ( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
//	bool bHasVolume = !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty();
	LayerSurface* layer = NULL;
	FSSurface* surf = NULL;
	if ( bHasSurface )
	{
		layer = ( LayerSurface* )( void* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() );
		if ( layer )
		{
			m_sliderOpacity->SetValue( (int)( layer->GetProperties()->GetOpacity() * 100 ) );
			double* rgb = layer->GetProperties()->GetColor();
			m_colorPicker->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
			rgb = layer->GetProperties()->GetEdgeColor();
			m_colorPickerEdge->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
			m_textFileName->ChangeValue( layer->GetFileName() );
			m_textFileName->SetInsertionPointEnd();
			m_textFileName->ShowPosition( m_textFileName->GetLastPosition() );
			m_spinEdgeThickness->SetValue( layer->GetProperties()->GetEdgeThickness() );
			surf = layer->GetSourceSurface();
		}
		
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
		lc->SetActiveLayer( ( Layer* )m_listBoxLayers->GetClientData( m_listBoxLayers->GetSelection() ) );
	}
	MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
	m_btnDelete->Enable( bHasSurface && !mainWnd->IsProcessing() );	
	m_btnSurfaceMain->Enable( bHasSurface );
	m_btnSurfaceInflated->Enable( bHasSurface && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceInflated ) );
	m_btnSurfaceWhite->Enable	( bHasSurface && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceWhite ) );
	m_btnSurfacePial->Enable	( bHasSurface && surf && surf->IsSurfaceLoaded( FSSurface::SurfacePial ) );
	m_btnSurfaceOriginal->Enable( bHasSurface && surf && surf->IsSurfaceLoaded( FSSurface::SurfaceOriginal ) );
	m_colorPicker->Enable( layer );
	m_colorPickerEdge->Enable( layer );
	m_choiceVector->Enable( layer );
	
	m_choiceVector->Clear();
	m_choiceVector->Append( "None" );
	if ( surf )
	{
		for ( int i = 0; i < surf->GetNumberOfVectorSets(); i++ )
		{
			m_choiceVector->Append( surf->GetVectorSetName( i ) );
		}
	}
	m_choiceVector->Append( "Load vector data..." );
	m_choiceVector->SetSelection( surf ? 1 + surf->GetActiveVector() : 0 );
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

void PanelSurface::UpdateLayerList( Layer* layer )
{
	LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
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

void PanelSurface::OnButtonDelete( wxCommandEvent& event )
{
	LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
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

void PanelSurface::OnButtonLoad( wxCommandEvent& event )
{
	MainWindow::GetMainWindowPointer()->LoadSurface();
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
		surf->GetProperties()->SetColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 ); 
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


void PanelSurface::OnButtonSurfaceMain( wxCommandEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->SetActiveSurface( FSSurface::SurfaceMain );
	}
}

void PanelSurface::OnButtonSurfaceInflated( wxCommandEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->SetActiveSurface( FSSurface::SurfaceInflated );
	}
}

void PanelSurface::OnButtonSurfaceWhite( wxCommandEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->SetActiveSurface( FSSurface::SurfaceWhite );
	}
}

void PanelSurface::OnButtonSurfacePial( wxCommandEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->SetActiveSurface( FSSurface::SurfacePial );
	}
}

void PanelSurface::OnButtonSurfaceOriginal( wxCommandEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->SetActiveSurface( FSSurface::SurfaceOriginal );
	}
}

void PanelSurface::OnSurfaceCloseUpdateUI( wxUpdateUIEvent& event )
{
	event.Enable( m_listBoxLayers->GetSelection() != wxNOT_FOUND );
}

void PanelSurface::OnSpinEdgeThickness( wxSpinEvent& event )
{
	LayerSurface* surf = ( LayerSurface* )MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
	if ( surf )
	{
		surf->GetProperties()->SetEdgeThickness( event.GetInt() );
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
