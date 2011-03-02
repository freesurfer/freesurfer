/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.17 $
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



#include "ToolWindowMeasure.h"
#include <wx/wx.h>
#include <wx/clrpicker.h>
#include <wx/config.h>
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include <wx/listbox.h>
#include <wx/clipbrd.h>
#include <wx/ffile.h>
#include "MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "BrushProperty.h"
#include "Interactor2DMeasure.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "Region2D.h"
#include "SurfaceRegion.h"
#include "SurfaceRegionGroups.h"

BEGIN_EVENT_TABLE( ToolWindowMeasure, wxFrame )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_LINE" ),      ToolWindowMeasure::OnActionMeasureLine )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_LINE" ),      ToolWindowMeasure::OnActionMeasureLineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_RECT" ),      ToolWindowMeasure::OnActionMeasureRectangle )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_RECT" ),      ToolWindowMeasure::OnActionMeasureRectangleUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_POLYLINE" ),  ToolWindowMeasure::OnActionMeasurePolyline )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_POLYLINE" ),  ToolWindowMeasure::OnActionMeasurePolylineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_SPLINE" ),    ToolWindowMeasure::OnActionMeasureSpline )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_SPLINE" ),    ToolWindowMeasure::OnActionMeasureSplineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_SURFACE" ),   ToolWindowMeasure::OnActionMeasureSurfaceRegion )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_SURFACE" ),   ToolWindowMeasure::OnActionMeasureSurfaceRegionUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_LABEL" ),     ToolWindowMeasure::OnActionMeasureLabel )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_LABEL" ),     ToolWindowMeasure::OnActionMeasureLabelUpdateUI )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_COPY" ),              ToolWindowMeasure::OnButtonCopy )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_EXPORT" ),            ToolWindowMeasure::OnButtonExport )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE" ),              ToolWindowMeasure::OnButtonSave )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE_ALL" ),          ToolWindowMeasure::OnButtonSaveAll )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_LOAD" ),              ToolWindowMeasure::OnButtonLoad )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_UPDATE" ),            ToolWindowMeasure::OnButtonUpdate )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_ID"),                   ToolWindowMeasure::OnSpinId )      
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_GROUP"),                   ToolWindowMeasure::OnSpinGroup )      
  EVT_COLOURPICKER_CHANGED ( XRCID( "ID_COLORPICKER_GROUP" ),     ToolWindowMeasure::OnColorGroup )
  
  EVT_SHOW      ( ToolWindowMeasure::OnShow )
  EVT_CLOSE     ( ToolWindowMeasure::OnClose )

END_EVENT_TABLE()


ToolWindowMeasure::ToolWindowMeasure( wxWindow* parent ) : Listener( "ToolWindowMeasure" )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_TOOLWINDOW_MEASURE") );
  m_toolbar     = XRCCTRL( *this, "ID_TOOLBAR_MEASURE", wxToolBar );
  m_textStats   = XRCCTRL( *this, "ID_TEXT_STATS",      wxTextCtrl );
  m_btnCopy     = XRCCTRL( *this, "ID_BUTTON_COPY",     wxButton );
  m_btnExport   = XRCCTRL( *this, "ID_BUTTON_EXPORT",   wxButton );
  m_btnSave     = XRCCTRL( *this, "ID_BUTTON_SAVE",     wxButton );
  m_btnSaveAll  = XRCCTRL( *this, "ID_BUTTON_SAVE_ALL", wxButton );
  m_btnLoad     = XRCCTRL( *this, "ID_BUTTON_LOAD",     wxButton );
  m_btnUpdate   = XRCCTRL( *this, "ID_BUTTON_UPDATE",   wxButton );
  m_spinId      = XRCCTRL( *this, "ID_SPIN_ID",         wxSpinCtrl );
  m_spinGroup   = XRCCTRL( *this, "ID_SPIN_GROUP",      wxSpinCtrl );
  m_colorPickerGroup = XRCCTRL( *this, "ID_COLORPICKER_GROUP",         wxColourPickerCtrl );
  
  m_widgets2D.push_back( m_btnCopy );
  m_widgets2D.push_back( m_btnExport );
  
  m_widgets3D.push_back( m_btnSave );
  m_widgets3D.push_back( m_btnSaveAll );
  m_widgets3D.push_back( m_btnLoad );
  m_widgets3D.push_back( m_spinId );
  m_widgets3D.push_back( m_spinGroup );
  m_widgets3D.push_back( m_colorPickerGroup );
  m_widgets3D.push_back( XRCCTRL( *this, "ID_STATIC_ID",    wxStaticText ) );
  m_widgets3D.push_back( XRCCTRL( *this, "ID_STATIC_GROUP",    wxStaticText ) );
  
  m_region = NULL;
  m_surfaceRegion = NULL;
  m_bToUpdateWidgets = true;
}

ToolWindowMeasure::~ToolWindowMeasure()
{}


void ToolWindowMeasure::OnClose( wxCloseEvent& event)
{
  Hide(); 
  MainWindow* mainwnd = MainWindow::GetMainWindowPointer();
  if ( mainwnd->IsShown() )
    mainwnd->SetMode( 0 );
}

void ToolWindowMeasure::OnShow( wxShowEvent& event )
{
//#if wxCHECK_VERSION(2,9,0)
#if wxVERSION_NUMBER > 2900  
  if ( event.IsShown() )
#else
  if ( event.GetShow() )
#endif  
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x = config->Read( _T("/ToolWindowMeasure/PosX"), 0L );
      int y = config->Read( _T("/ToolWindowMeasure/PosY"), 0L );
      if ( x == 0 && y == 0 )
        Center();
      else
        Move( x, y );
    }
    
 //   SetClientSize( GetClientSize().GetWidth(), m_toolbar->GetSize().GetHeight() );
  }
  else
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x, y;
      GetPosition( &x, &y );
      config->Write( _T("/ToolWindowMeasure/PosX"), (long) x );
      config->Write( _T("/ToolWindowMeasure/PosY"), (long) y );
    }
  }
}

void ToolWindowMeasure::SetRegion( Region2D* reg )
{
  m_region = reg;
  if ( m_region )
  {
    m_region->AddListener( this );
    m_surfaceRegion = NULL;
    RenderView* view = MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
    if ( view->GetAction() == Interactor::MM_SurfaceRegion )
      MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Line );
  }  
  UpdateWidgets();
}

void ToolWindowMeasure::SetSurfaceRegion( SurfaceRegion* reg )
{
  m_surfaceRegion = reg;
  if ( m_surfaceRegion )
  {
    m_surfaceRegion->AddListener( this );
    m_region = NULL;
    MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_SurfaceRegion );
  }
  UpdateWidgets();
}

void ToolWindowMeasure::UpdateWidgets( )
{
  m_bToUpdateWidgets = true;
}

wxString ToolWindowMeasure::GetLabelStats()
{
  wxString strg;
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerMRI* label = NULL, *mri = NULL;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)lc->GetLayer( i ) )->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
    {
      label = ( (LayerMRI*)lc->GetLayer( i ) );
      break;
    }
  }
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)lc->GetLayer( i ) )->GetProperties()->GetColorMap() != LayerPropertiesMRI::LUT )
    {
      mri = ( (LayerMRI*)lc->GetLayer( i ) );
      break;
    }
  }
  if ( label && mri )
  {
    int nPlane = MainWindow::GetMainWindowPointer()->GetMainViewId();
    std::vector<int> ids, numbers;
    std::vector<double> means, sds;
    mri->GetLabelStats( label, nPlane, ids, numbers, means, sds );
    strg << "Id \tCount \tMean \t+/-SD\n";
    for ( size_t i = 0; i < ids.size(); i++ )
    {
      wxString snum, smean;
      snum << numbers[i];
      smean << means[i];
      if ( snum.size() < 4 )
        snum.Pad( (4-snum.size())*2 );
      if ( smean.size() < 4 )
        smean.Pad( (4-smean.size())*2 );
      strg << ids[i] << " \t" << snum << " \t" << smean << " \t" << sds[i] << "\n";
    }
  }
  
  return strg;
}

void ToolWindowMeasure::DoUpdateWidgets()
{
  wxString strg;
  RenderView2D* view = (RenderView2D*)MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  if ( view->GetAction() == Interactor::MM_Label )
  {
    strg = GetLabelStats();
  }
  else if ( m_region )
  {
    wxArrayString strgs = m_region->GetLongStats();
    for ( size_t i = 0; i < strgs.size(); i++ )
      strg += strgs[i] + "\n";   
  }
  m_textStats->ChangeValue( strg ); 
  m_btnCopy->Enable( !strg.IsEmpty() );
  m_btnExport->Enable( !strg.IsEmpty() );
  m_btnUpdate->Show( view->GetAction() == Interactor::MM_Label );
  
  for ( size_t i = 0; i < m_widgets3D.size(); i++ )
    m_widgets3D[i]->Show( view->GetAction() == Interactor::MM_SurfaceRegion );
  
  for ( size_t i = 0; i < m_widgets2D.size(); i++ )
    m_widgets2D[i]->Show( view->GetAction() != Interactor::MM_SurfaceRegion );
  
  if ( m_surfaceRegion )
  {
    m_spinId->SetValue( m_surfaceRegion->GetId() );
    m_spinGroup->SetValue( m_surfaceRegion->GetGroup() );
    m_spinGroup->SetRange( 1, 
                           m_surfaceRegion->GetMRI()->GetSurfaceRegionGroups()
                               ->GetGroupIdRange( m_surfaceRegion ) );
    m_colorPickerGroup->SetColour( m_surfaceRegion->GetColor() );
  }
  
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  bool bSurfaceRegionValid = ( mri && mri->GetProperties()->GetShowAsContour() && mri->GetNumberOfSurfaceRegions() > 0 );
  if ( bSurfaceRegionValid )
    m_spinId->SetRange( 1, mri->GetNumberOfSurfaceRegions() );
  
  m_btnSave->Enable( m_surfaceRegion && bSurfaceRegionValid );
  m_spinId->Enable( m_surfaceRegion && bSurfaceRegionValid );
  m_spinGroup->Enable( m_surfaceRegion && bSurfaceRegionValid );
  m_colorPickerGroup->Enable( m_surfaceRegion && bSurfaceRegionValid );
  m_btnSaveAll->Enable( bSurfaceRegionValid );
  
  XRCCTRL( *this, "ID_PANEL_HOLDER",  wxPanel )->Layout();
  Layout();
  m_bToUpdateWidgets = false;
}

void ToolWindowMeasure::OnInternalIdle()
{
  wxFrame::OnInternalIdle();
  
  if ( m_bToUpdateWidgets )
    DoUpdateWidgets();
}

void ToolWindowMeasure::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "RegionStatsUpdated" )
  {
    if ( m_region == iData || m_region == sender )
      UpdateWidgets();
  }
  else if ( iMsg == "RegionSelected" )
  {
    SetRegion( (Region2D*)iData );
  }
  else if ( iMsg == "RegionRemoved" )
  {
    if ( m_region == iData )
      SetRegion( NULL );
  }
  else if ( iMsg == "SurfaceRegionSelected" )
  {
    SetSurfaceRegion( (SurfaceRegion*)iData );
  }
  else if ( iMsg == "SurfaceRegionRemoved" )
  {
    if ( m_surfaceRegion == iData )
      SetSurfaceRegion( NULL );
  }
  else if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || iMsg == "LayerEdited" || iMsg == "SurfaceRegionColorChanged" )
  {
    UpdateWidgets();
  }
}

void ToolWindowMeasure::OnActionMeasureLine( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Line );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasureLineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_Line );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowMeasure::OnActionMeasureRectangle( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Rectangle );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasureRectangleUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_Rectangle );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowMeasure::OnActionMeasurePolyline( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Polyline );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasurePolylineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_Polyline );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowMeasure::OnActionMeasureSpline( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Spline );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasureSplineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_Spline );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowMeasure::OnActionMeasureSurfaceRegion( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_SurfaceRegion );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasureSurfaceRegionUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_SurfaceRegion );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowMeasure::OnActionMeasureLabel( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor::MM_Label );
  UpdateWidgets();
}

void ToolWindowMeasure::OnActionMeasureLabelUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor::MM_Label );
  LayerCollection* col = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  bool bLabelExist = false;
  for ( int i = 0; i < col->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)col->GetLayer( i ) )->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
    {
      bLabelExist = true;
      break;
    }
  }
  event.Enable( MainWindow::GetMainWindowPointer()->GetMainViewId() < 3 &&
      view->GetInteractionMode() == RenderView2D::IM_Measure &&
      col->GetNumberOfLayers() > 1 && bLabelExist );
}

void ToolWindowMeasure::OnButtonCopy( wxCommandEvent& event )
{
  wxString output = m_textStats->GetValue();
  if (wxTheClipboard->Open())
  {
    wxTheClipboard->SetData( new wxTextDataObject( output ) );
    wxTheClipboard->Close();
  }
}

void ToolWindowMeasure::OnButtonExport( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Export to file"), _(""), _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxString fn = dlg.GetPath();
    if ( !m_textStats->SaveFile( fn ) )
    {
      wxMessageDialog msg_dlg( this, wxString("Can not write to file ") + fn, 
                           _("Error"), wxOK );
      msg_dlg.ShowModal();
    } 
  }
}

void ToolWindowMeasure::OnSpinId( wxSpinEvent& event )
{
  RenderView3D* view = ( RenderView3D* )MainWindow::GetMainWindowPointer()->GetRenderView( 3 );
  view->PickSelectRegion( event.GetInt() );
  UpdateWidgets();
}

void ToolWindowMeasure::OnButtonSave( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Save region to file"), _(""), _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( m_surfaceRegion && dlg.ShowModal() == wxID_OK )
  {
    wxString fn = dlg.GetPath();
    if ( !m_surfaceRegion->Write( fn ) )
    {
      wxMessageDialog msg_dlg( this, wxString("Can not write to file ") + fn, 
                               _("Error"), wxOK );
      msg_dlg.ShowModal();
    } 
  }
}

void ToolWindowMeasure::OnButtonSaveAll( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Save all regions to file"), _(""), _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxString fn = dlg.GetPath();
    LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
    if ( mri && !mri->SaveAllSurfaceRegions( fn ) )
    {
      wxMessageDialog msg_dlg( this, wxString("Can not write to file ") + fn, 
                               _("Error"), wxOK );
      msg_dlg.ShowModal();
    } 
  }
}


void ToolWindowMeasure::OnButtonLoad( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Load region(s) from file"), _(""), _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_OPEN );
  
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( mri && dlg.ShowModal() == wxID_OK )
  {
    wxString fn = dlg.GetPath();
    if ( !mri->LoadSurfaceRegions( fn ) )
    {
      wxMessageDialog msg_dlg( this, wxString("Can not load file ") + fn, 
                               _("Error"), wxOK );
      msg_dlg.ShowModal();
    } 
    UpdateWidgets();
  }
}

void ToolWindowMeasure::OnButtonUpdate( wxCommandEvent& event )
{
  UpdateWidgets();
}

void ToolWindowMeasure::OnSpinGroup( wxSpinEvent& event )
{
  if ( m_surfaceRegion )
  {
    m_surfaceRegion->SetGroup( event.GetInt() );
  }
}
  
void ToolWindowMeasure::OnColorGroup( wxColourPickerEvent& event )
{
  if ( m_surfaceRegion )
  {
    m_surfaceRegion->GetMRI()->GetSurfaceRegionGroups()
        ->SetGroupColor( m_surfaceRegion->GetGroup(), event.GetColour() );
  }
}
