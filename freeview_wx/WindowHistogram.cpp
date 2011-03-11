/**
 * @file  WindowHistogram.h
 * @brief Main window.
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

#include <iostream>
#include <wx/wx.h>
#include "WindowHistogram.h"
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
#include <vtkImageData.h>

BEGIN_EVENT_TABLE( WindowHistogram, wxFrame )
	EVT_CLOSE		  ( WindowHistogram::OnClose )
	EVT_BUTTON		( XRCID( "ID_BUTTON_CLOSE" ), 	WindowHistogram::OnButtonClose )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_LAYER" ),   WindowHistogram::OnChoiceLayer )
END_EVENT_TABLE()

WindowHistogram::WindowHistogram( wxWindow* parent ) : Listener( "WindowHistogram" )
{
	wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_FRAME_HISTOGRAM") );
	m_histogramWidget = XRCCTRL( *this, "ID_HISTOGRAM_WIDGET", wxHistogramWidget );
	m_choiceLayer = XRCCTRL( *this, "ID_CHOICE_LAYER", wxChoice );
 
	wxConfigBase* config = wxConfigBase::Get();
	if ( config )
	{
		int x = config->Read( _T("/HistogramWindow/PosX"), 280L );
		int y = config->Read( _T("/HistogramWindow/PosY"), 30L );
		int w = config->Read( _T("/HistogramWindow/Width"), 380L );
		int h = config->Read( _T("/HistogramWindow/Height"), 520L );
		SetSize( w, h );
		int x1 = 0, y1 = 0;
		if ( parent )
			parent->GetPosition( &x1, &y1 );
		Move( x1 + x, y1 + y );
        
    m_histogramWidget->SetNumberOfBins( config->Read( _T("/HistogramWindow/NumberOfBins"), 100L ) );
	}
	
	m_activeLayer = NULL;
	
	MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->AddListener( this );
	MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->AddListener( this );
	
	UpdateUI();
}

void WindowHistogram::OnClose( wxCloseEvent& event )
{
	wxConfigBase* config = wxConfigBase::Get();
	if ( config && !IsIconized() )
	{
		int x, x2, y, y2, w, h;
		GetParent()->GetPosition( &x2, &y2 );
		GetPosition( &x, &y );
		GetSize( &w, &h );
		config->Write( _T("/HistogramWindow/PosX"), (long) x - x2 );
		config->Write( _T("/HistogramWindow/PosY"), (long) y - y2 );
		config->Write( _T("/HistogramWindow/Width"), (long) w );
		config->Write( _T("/HistogramWindow/Height"), (long) h );
	}

	Hide();
}


void WindowHistogram::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
	if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" )
	{
	//	UpdateUI();
	}
}

void WindowHistogram::ShowWindow( Layer* default_layer )
{
	if ( m_activeLayer != default_layer )
	{
		UpdateUI();
		UpdateHistogram();
	}
	
	Show();
}

// unfinished!
void WindowHistogram::UpdateUI()
{
	m_choiceLayer->Clear();
	LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
	LayerCollection* lc_surf = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
	std::vector<Layer*> layers = lc_mri->GetLayers();
	for ( int i = 0; i < lc_surf->GetNumberOfLayers(); i++ )
	{
		layers.push_back( lc_surf->GetLayer( i ) );
	}
	int n = -1;
	for ( size_t i = 0; i < layers.size(); i++ )
	{
    m_choiceLayer->Append( wxString::FromAscii( layers[i]->GetName() ), ( void* )layers[i] );
		if ( layers[i] == m_activeLayer )
			n = i;
	}
	
	if ( n < 0 && layers.size() > 0 )
	{
		m_activeLayer = layers[0];
        UpdateHistogram();
		n = 0;
	}
	
	if ( n >= 0 )		// at least has a layer
	{
		m_choiceLayer->SetSelection( n );
	}
	else
	{
		Close();
	}
}

void WindowHistogram::UpdateHistogram()
{
	if ( !m_activeLayer )
		return;
	
	if ( m_activeLayer->IsTypeOf( "MRI" ) )
	{
    LayerMRI* mri = (LayerMRI*)m_activeLayer;
    vtkImageData* image = mri->GetImageData();
    int nSize = image->GetNumberOfPoints();
    switch ( image->GetScalarType() )
    { 
    case VTK_INT:
      m_histogramWidget->SetInputData( ( int* )image->GetScalarPointer(), nSize );
      break;
    case VTK_FLOAT:
      m_histogramWidget->SetInputData( ( float* )image->GetScalarPointer(), nSize );
      break; 
    case VTK_SHORT:
      m_histogramWidget->SetInputData( ( short* )image->GetScalarPointer(), nSize );
      break;
    case VTK_UNSIGNED_SHORT:
      m_histogramWidget->SetInputData( ( unsigned short* )image->GetScalarPointer(), nSize );
      break;
    case VTK_DOUBLE:
      m_histogramWidget->SetInputData( ( double* )image->GetScalarPointer(), nSize );
      break;
    case VTK_UNSIGNED_CHAR:
      m_histogramWidget->SetInputData( ( unsigned char* )image->GetScalarPointer(), nSize );
      break;
    case VTK_LONG:
      m_histogramWidget->SetInputData( ( long* )image->GetScalarPointer(), nSize );
      break;
    default:
      std::cerr << "Unsupported data type" << std::endl;
      break;
    }      
	}
}

void WindowHistogram::OnChoiceLayer( wxCommandEvent& event )
{
  m_activeLayer = (Layer*)(void*)m_choiceLayer->GetClientData( event.GetSelection() );
  UpdateUI();
  UpdateHistogram();
}
