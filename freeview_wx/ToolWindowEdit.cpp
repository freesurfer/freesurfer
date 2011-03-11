/**
 * @file  ToolWindowEdit.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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


#include "ToolWindowEdit.h"
#include <wx/wx.h>
#include <wx/config.h>
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include <wx/clrpicker.h>
#include "MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "BrushProperty.h"
#include "Interactor2DROIEdit.h"
#include "Interactor2DVoxelEdit.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LayerROI.h"
#include "LayerDTI.h"
#include "Contour2D.h"

BEGIN_EVENT_TABLE( ToolWindowEdit, wxFrame )
// EVT_BUTTON   ( wxID_OK,          ToolWindowEdit::OnOK )
// EVT_BUTTON   ( XRCID( wxT( "ID_BUTTON_VECTOR_FILE" ) ), ToolWindowEdit::OnButtonVector )
// EVT_BUTTON   ( XRCID( wxT( "ID_BUTTON_FA_FILE" ) ),  ToolWindowEdit::OnButtonFA )
// EVT_COMBOBOX  ( XRCID( wxT( "ID_COMBO_FA_FILE" ) ),   ToolWindowEdit::OnComboFASelectionChanged )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_FREEHAND" ),  ToolWindowEdit::OnActionVoxelFreehand )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_FREEHAND" ),  ToolWindowEdit::OnActionVoxelFreehandUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_FILL" ),      ToolWindowEdit::OnActionVoxelFill )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_FILL" ),      ToolWindowEdit::OnActionVoxelFillUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_POLYLINE" ),  ToolWindowEdit::OnActionVoxelPolyline )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_POLYLINE" ),  ToolWindowEdit::OnActionVoxelPolylineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_LIVEWIRE" ),  ToolWindowEdit::OnActionVoxelLivewire )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_LIVEWIRE" ),  ToolWindowEdit::OnActionVoxelLivewireUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_EYEDROP" ),   ToolWindowEdit::OnActionVoxelColorPicker )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_EYEDROP" ),   ToolWindowEdit::OnActionVoxelColorPickerUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_CONTOUR" ),   ToolWindowEdit::OnActionVoxelContour )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_CONTOUR" ),   ToolWindowEdit::OnActionVoxelContourUpdateUI )
  
  EVT_MENU      ( XRCID( "ID_ACTION_VOXEL_LIVEWIRE" ),  ToolWindowEdit::OnActionVoxelLivewire )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_VOXEL_LIVEWIRE" ),  ToolWindowEdit::OnActionVoxelLivewireUpdateUI )
  
  EVT_MENU      ( XRCID( "ID_ACTION_ROI_FREEHAND" ),    ToolWindowEdit::OnActionROIFreehand )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_ROI_FREEHAND" ),    ToolWindowEdit::OnActionROIFreehandUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_ROI_FILL" ),        ToolWindowEdit::OnActionROIFill )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_ROI_FILL" ),        ToolWindowEdit::OnActionROIFillUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_ROI_POLYLINE" ),    ToolWindowEdit::OnActionROIPolyline )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_ROI_POLYLINE" ),    ToolWindowEdit::OnActionROIPolylineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_ROI_LIVEWIRE" ),    ToolWindowEdit::OnActionROILivewire )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_ROI_LIVEWIRE" ),    ToolWindowEdit::OnActionROILivewireUpdateUI )
  
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_SIZE" ),        ToolWindowEdit::OnSpinBrushSize )
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_BRUSH_TOLERANCE" ),   ToolWindowEdit::OnSpinBrushTolerance )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_REFERENCE" ),       ToolWindowEdit::OnChoiceBrushTemplate )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_DRAW_CONNECTED" ),   ToolWindowEdit::OnCheckDrawConnectedOnly )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_DRAW_RANGE" ),       ToolWindowEdit::OnCheckDrawRange )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_EXCLUDE_RANGE" ),    ToolWindowEdit::OnCheckExcludeRange )
  EVT_TEXT      ( XRCID( "ID_EDIT_DRAW_RANGE_LOW" ),    ToolWindowEdit::OnEditDrawRangeLow )
  EVT_TEXT      ( XRCID( "ID_EDIT_DRAW_RANGE_HIGH" ),   ToolWindowEdit::OnEditDrawRangeHigh )
  EVT_TEXT      ( XRCID( "ID_EDIT_EXCLUDE_RANGE_LOW" ), ToolWindowEdit::OnEditExcludeRangeLow )
  EVT_TEXT      ( XRCID( "ID_EDIT_EXCLUDE_RANGE_HIGH" ),  ToolWindowEdit::OnEditExcludeRangeHigh )
  
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_SMOOTH" ),           ToolWindowEdit::OnCheckSmooth )
  EVT_TEXT      ( XRCID( "ID_EDIT_SMOOTH_SD" ),         ToolWindowEdit::OnEditSmoothSD )
  EVT_TEXT      ( XRCID( "ID_EDIT_CONTOUR_VALUE" ),     ToolWindowEdit::OnEditContourValue )
  EVT_COLOURPICKER_CHANGED  ( XRCID( "ID_COLORPICKER_CONTOUR" ),  ToolWindowEdit::OnColorContour )

  EVT_SHOW      ( ToolWindowEdit::OnShow )
  EVT_CLOSE     ( ToolWindowEdit::OnClose )

END_EVENT_TABLE()


ToolWindowEdit::ToolWindowEdit( wxWindow* parent ) : Listener( "ToolWindowMeasure" ),
    m_editSmoothSD( NULL ),
    m_bToUpdateTools( false )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_TOOLWINDOW_EDIT") );
  m_toolbarVoxelEdit    = XRCCTRL( *this, "ID_TOOLBAR_VOXEL_EDIT",    wxToolBar );
  m_toolbarROIEdit      = XRCCTRL( *this, "ID_TOOLBAR_ROI_EDIT",      wxToolBar );
  m_spinBrushSize       = XRCCTRL( *this, "ID_SPIN_BRUSH_SIZE",       wxSpinCtrl );
  m_spinBrushTolerance  = XRCCTRL( *this, "ID_SPIN_BRUSH_TOLERANCE",  wxSpinCtrl );
  m_choiceTemplate      = XRCCTRL( *this, "ID_CHOICE_REFERENCE",      wxChoice );
  m_checkDrawRange      = XRCCTRL( *this, "ID_CHECK_DRAW_RANGE",      wxCheckBox );
  m_checkExcludeRange   = XRCCTRL( *this, "ID_CHECK_EXCLUDE_RANGE",   wxCheckBox );
  m_editDrawRangeLow    = XRCCTRL( *this, "ID_EDIT_DRAW_RANGE_LOW",   wxTextCtrl );
  m_editDrawRangeHigh   = XRCCTRL( *this, "ID_EDIT_DRAW_RANGE_HIGH",  wxTextCtrl );
  m_editExcludeRangeLow     = XRCCTRL( *this, "ID_EDIT_EXCLUDE_RANGE_LOW",  wxTextCtrl );
  m_editExcludeRangeHigh    = XRCCTRL( *this, "ID_EDIT_EXCLUDE_RANGE_HIGH", wxTextCtrl );
  m_checkDrawConnectedOnly  = XRCCTRL( *this, "ID_CHECK_DRAW_CONNECTED",    wxCheckBox );
  m_checkSmooth         = XRCCTRL( *this, "ID_CHECK_SMOOTH",          wxCheckBox );
  m_editSmoothSD        = XRCCTRL( *this, "ID_EDIT_SMOOTH_SD",        wxTextCtrl );
  m_editContourValue    = XRCCTRL( *this, "ID_EDIT_CONTOUR_VALUE",    wxTextCtrl );
  m_colorPickerContour  = XRCCTRL( *this, "ID_COLORPICKER_CONTOUR",   wxColourPickerCtrl );
  
  m_widgetsBrushSize.push_back( XRCCTRL( *this, "ID_STATIC_BRUSH_SIZE", wxStaticText ) );
  m_widgetsBrushSize.push_back( m_spinBrushSize );
  
  m_widgetsReference.push_back( XRCCTRL( *this, "ID_STATIC_REFERENCE", wxStaticText ) );
  m_widgetsReference.push_back( m_choiceTemplate );
  
  m_widgetsTolerance.push_back( m_spinBrushTolerance );
  m_widgetsTolerance.push_back( XRCCTRL( *this, "ID_STATIC_TOLERANCE", wxStaticText ) );
  m_widgetsTolerance.push_back( XRCCTRL( *this, "ID_STATIC_PERCENTAGE", wxStaticText ) );
  
  m_widgetsConstrain.push_back( m_checkDrawConnectedOnly );  
  m_widgetsConstrain.push_back( m_checkDrawRange );
  m_widgetsConstrain.push_back( m_editDrawRangeLow );
  m_widgetsConstrain.push_back( m_editDrawRangeHigh );
  m_widgetsConstrain.push_back( m_checkExcludeRange );
  m_widgetsConstrain.push_back( m_editExcludeRangeLow );
  m_widgetsConstrain.push_back( m_editExcludeRangeHigh );
  m_widgetsConstrain.push_back( XRCCTRL( *this, "ID_STATIC_DRAW_LOW", wxStaticText ) );
  m_widgetsConstrain.push_back( XRCCTRL( *this, "ID_STATIC_DRAW_HIGH", wxStaticText ) );
  m_widgetsConstrain.push_back( XRCCTRL( *this, "ID_STATIC_EXCLUDE_LOW", wxStaticText ) );
  m_widgetsConstrain.push_back( XRCCTRL( *this, "ID_STATIC_EXCLUDE_HIGH", wxStaticText ) );
  
  m_widgetsSmooth.push_back( m_checkSmooth );
  m_widgetsSmooth.push_back( m_editSmoothSD );
  m_widgetsSmooth.push_back( XRCCTRL( *this, "ID_STATIC_SMOOTH_SD", wxStaticText ) );
  
  m_widgetsContour.push_back( XRCCTRL( *this, "ID_STATIC_NOTES_CONTOUR", wxStaticText ) );
  m_widgetsContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_VALUE", wxStaticText ) );
  m_widgetsContour.push_back( m_editContourValue );
  m_widgetsContour.push_back( XRCCTRL( *this, "ID_STATIC_CONTOUR_COLOR", wxStaticText ) );
  m_widgetsContour.push_back( m_colorPickerContour );
}

ToolWindowEdit::~ToolWindowEdit()
{}


void ToolWindowEdit::ShowWidgets( std::vector<wxWindow*>& list, bool bShow )
{
  for ( size_t i = 0; i < list.size(); i++ )
  {
    list[i]->Show( bShow );
  }
}

void ToolWindowEdit::OnClose( wxCloseEvent& event )
{
  Hide();
  MainWindow* mainwnd = MainWindow::GetMainWindowPointer();
  if ( mainwnd->IsShown() )
    mainwnd->SetMode( 0 );
}

void ToolWindowEdit::OnShow( wxShowEvent& event )
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
      int x = config->Read( _T("/ToolWindowEdit/PosX"), 0L );
      int y = config->Read( _T("/ToolWindowEdit/PosY"), 0L );
      if ( x == 0 && y == 0 )
        Center();
      else
        Move( x, y );
    }
    for ( int i = 0; i < 3; i++ )
    {
      RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( i );
      view->GetContour2D()->AddListener( this );
    }
  }
  else
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x, y;
      GetPosition( &x, &y );
      config->Write( _T("/ToolWindowEdit/PosX"), (long) x );
      config->Write( _T("/ToolWindowEdit/PosY"), (long) y );
    } 
  }
  MainWindow::GetMainWindowPointer()->SetFocus();
}

void ToolWindowEdit::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "ContourValueChanged" )
  {
    UpdateTools();
  }
}

void ToolWindowEdit::UpdateTools()
{
  m_bToUpdateTools = true;
}

void ToolWindowEdit::DoUpdateTools()
{
  int nViewId = MainWindow::GetMainWindowPointer()->GetActiveViewId();
  if ( nViewId < 0 || nViewId > 2 )
    nViewId = 0;
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( nViewId );
  
  bool bVoxelEditVisible = m_toolbarVoxelEdit->IsShown();
  bool bROIEditVisible = m_toolbarROIEdit->IsShown();
  if ( bVoxelEditVisible != (view->GetInteractionMode() == RenderView2D::IM_VoxelEdit) ||
       bROIEditVisible != (view->GetInteractionMode() == RenderView2D::IM_ROIEdit) )
  {
    m_toolbarVoxelEdit ->Show( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit );
    m_toolbarROIEdit ->Show( view->GetInteractionMode() == RenderView2D::IM_ROIEdit );

    XRCCTRL( *this, "ID_PANEL_HOLDER", wxPanel )->Layout();
  }

// XRCCTRL( *m_toolbarBrush, "ID_STATIC_BRUSH_SIZE", wxStaticText )->Enable( m_viewAxial->GetAction() != Interactor2DROIEdit::EM_Fill );
  m_spinBrushSize->Enable( view->GetAction() != Interactor2DROIEdit::EM_Fill );
  m_spinBrushTolerance->Enable( view->GetAction() == Interactor2DROIEdit::EM_Fill );
// choiceTemplate->Enable( checkTemplate->IsChecked() && m_viewAxial->GetAction() == Interactor2DROIEdit::EM_Fill );
// XRCCTRL( *m_toolbarBrush, "ID_STATIC_BRUSH_TOLERANCE", wxStaticText )->Enable( checkTemplate->IsChecked() ); //&& m_viewAxial->GetAction() == Interactor2DROIEdit::EM_Fill );
// XRCCTRL( *m_toolbarBrush, "ID_SPIN_BRUSH_TOLERANCE", wxSpinCtrl )->Enable( checkTemplate->IsChecked() );//&& m_viewAxial->GetAction() == Interactor2DROIEdit::EM_Fill );

  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  LayerVolumeBase* layer = bp->GetReferenceLayer();
// if ( m_choiceTemplate->GetSelection() != wxNOT_FOUND )
//  layer = ( LayerEditable* )(void*)m_choiceTemplate->GetClientData( m_choiceTemplate->GetSelection() );

  m_choiceTemplate->Clear();
  m_choiceTemplate->Append( _("None"), (void*)NULL );
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  int nSel = 0;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri = (LayerMRI*)lc->GetLayer( i );
    if ( layer == mri )
    {
      nSel = i+1;
    }

    m_choiceTemplate->Append( wxString::FromAscii( mri->GetName() ), (void*)mri );
  }
// if ( !lc->IsEmpty() )
  m_choiceTemplate->SetSelection( nSel );

  m_spinBrushSize->SetValue( bp->GetBrushSize() );
  m_spinBrushTolerance->SetValue( bp->GetBrushTolerance( ) );

  m_checkDrawConnectedOnly->SetValue( bp->GetDrawConnectedOnly() );
  m_checkDrawRange  ->SetValue( bp->GetDrawRangeEnabled() );
  m_checkExcludeRange  ->SetValue( bp->GetExcludeRangeEnabled() );

  m_editDrawRangeLow  ->Enable( bp->GetDrawRangeEnabled() );
  m_editDrawRangeHigh  ->Enable( bp->GetDrawRangeEnabled() );
  m_editExcludeRangeLow ->Enable( bp->GetExcludeRangeEnabled() );
  m_editExcludeRangeHigh ->Enable( bp->GetExcludeRangeEnabled() );

  double* range = bp->GetDrawRange();
  UpdateTextValue( m_editDrawRangeLow, range[0] );
  UpdateTextValue( m_editDrawRangeHigh, range[1] );
  range = bp->GetExcludeRange();
  UpdateTextValue( m_editExcludeRangeLow, range[0] );
  UpdateTextValue( m_editExcludeRangeHigh, range[1] );

  Contour2D* c2d = view->GetContour2D();
  m_checkSmooth->SetValue( c2d->GetSmooth() );
  UpdateTextValue( m_editSmoothSD, c2d->GetSmoothSD() );
  m_editSmoothSD->Enable( c2d->GetSmooth() );
  UpdateTextValue( m_editContourValue, c2d->GetContourValue() );
  
  double* rgb = c2d->GetContourColor();
  m_colorPickerContour->SetColour( wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
  
  int nAction = view->GetAction(); 
  ShowWidgets( m_widgetsBrushSize, nAction != Interactor2DVoxelEdit::EM_Contour &&
                                   nAction != Interactor2DVoxelEdit::EM_ColorPicker && 
                                   nAction != Interactor2DVoxelEdit::EM_Fill );
  ShowWidgets( m_widgetsReference, nAction == Interactor2DVoxelEdit::EM_Fill || 
                                   nAction == Interactor2DVoxelEdit::EM_Contour );  
  ShowWidgets( m_widgetsTolerance, nAction == Interactor2DVoxelEdit::EM_Fill );
  ShowWidgets( m_widgetsConstrain, nAction != Interactor2DVoxelEdit::EM_ColorPicker && 
      nAction != Interactor2DVoxelEdit::EM_Contour ); 
  ShowWidgets( m_widgetsSmooth, nAction == Interactor2DVoxelEdit::EM_Contour );
  ShowWidgets( m_widgetsContour, nAction == Interactor2DVoxelEdit::EM_Contour );
  
  m_bToUpdateTools = false;
  wxPanel* panel = XRCCTRL( *this, "ID_PANEL_HOLDER", wxPanel );
  panel->Layout();
  panel->Fit();
  Fit();
  Layout();
  MainWindow::GetMainWindowPointer()->NeedRedraw( 1 );
}

void ToolWindowEdit::UpdateTextValue( wxTextCtrl* ctrl, double dvalue )
{
  double dtemp;
  if ( ctrl->GetValue().ToDouble( &dtemp ) && dtemp == dvalue )
    return;
  
  wxString value_strg = ( wxString() << dvalue );
  if ( value_strg != ctrl->GetValue() && (value_strg + _(".")) != ctrl->GetValue() )
    ctrl->ChangeValue( value_strg );
}

void ToolWindowEdit::OnInternalIdle()
{
  wxFrame::OnInternalIdle();

  if ( m_bToUpdateTools )
    DoUpdateTools();
}

void ToolWindowEdit::OnActionVoxelFreehand( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_Freehand );
  MainWindow::GetMainWindowPointer()->GetBrushProperty()->SetReferenceLayer( NULL );
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelFreehandUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
               && view->GetAction() == Interactor2DVoxelEdit::EM_Freehand );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowEdit::OnActionVoxelPolyline( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_Polyline );
  MainWindow::GetMainWindowPointer()->GetBrushProperty()->SetReferenceLayer( NULL );
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelPolylineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
               && view->GetAction() == Interactor2DVoxelEdit::EM_Polyline );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowEdit::OnActionVoxelLivewire( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_Livewire );
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelLivewireUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
               && view->GetAction() == Interactor2DVoxelEdit::EM_Livewire );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowEdit::OnActionVoxelFill( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_Fill );
  MainWindow::GetMainWindowPointer()->GetBrushProperty()->SetReferenceLayer( NULL );
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelFillUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
               && view->GetAction() == Interactor2DVoxelEdit::EM_Fill );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowEdit::OnActionVoxelColorPicker( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_ColorPicker );
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelColorPickerUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
      && view->GetAction() == Interactor2DVoxelEdit::EM_ColorPicker );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowEdit::OnActionVoxelContour( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DVoxelEdit::EM_Contour );
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();  
  LayerMRI* layer = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( layer && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
    layer->GetProperties()->SetShowLabelOutline( true );
  
  if ( !bp->GetReferenceLayer() )
  {
    LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
    for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
    {
      LayerMRI* mri = (LayerMRI*)lc_mri->GetLayer( i );
      if ( mri->GetProperties()->GetColorMap() != LayerPropertiesMRI::LUT && mri->IsVisible() && mri != layer )
      {
        bp->SetReferenceLayer( mri );
        break;
      }
    }
  }
  UpdateTools();
}

void ToolWindowEdit::OnActionVoxelContourUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
      && view->GetAction() == Interactor2DVoxelEdit::EM_Contour );

  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer();
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_VoxelEdit
      && mri /*&& mri->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT*/ );
}

void ToolWindowEdit::OnActionROIFreehand( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DROIEdit::EM_Freehand );
  UpdateTools();
}

void ToolWindowEdit::OnActionROIFreehandUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
               && view->GetAction() == Interactor2DROIEdit::EM_Freehand );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowEdit::OnActionROIPolyline( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DROIEdit::EM_Polyline );
  UpdateTools();
}

void ToolWindowEdit::OnActionROIPolylineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
               && view->GetAction() == Interactor2DROIEdit::EM_Polyline );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowEdit::OnActionROILivewire( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DROIEdit::EM_Livewire );
  UpdateTools();
}

void ToolWindowEdit::OnActionROILivewireUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
               && view->GetAction() == Interactor2DROIEdit::EM_Livewire );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowEdit::OnActionROIFill( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DROIEdit::EM_Fill );
  UpdateTools();
}

void ToolWindowEdit::OnActionROIFillUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
               && view->GetAction() == Interactor2DROIEdit::EM_Fill );

  event.Enable( view->GetInteractionMode() == RenderView2D::IM_ROIEdit
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowEdit::OnSpinBrushSize( wxSpinEvent& event )
{
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetBrushSize( event.GetInt() );
  UpdateTools();
}

void ToolWindowEdit::OnSpinBrushTolerance( wxSpinEvent& event )
{
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetBrushTolerance( event.GetInt() );
  UpdateTools();
}

void ToolWindowEdit::OnChoiceBrushTemplate( wxCommandEvent& event )
{
  LayerVolumeBase* layer = (LayerVolumeBase*)(void*)m_choiceTemplate->GetClientData( event.GetSelection() );
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetReferenceLayer( layer );
  UpdateTools();
}

void ToolWindowEdit::OnCheckDrawConnectedOnly( wxCommandEvent& event )
{
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetDrawConnectedOnly( event.IsChecked() );

  UpdateTools();
}

void ToolWindowEdit::OnCheckDrawRange( wxCommandEvent& event )
{
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetDrawRangeEnabled( event.IsChecked() );

  UpdateTools();
}

void ToolWindowEdit::OnCheckExcludeRange( wxCommandEvent& event )
{
  BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
  bp->SetExcludeRangeEnabled( event.IsChecked() );
  UpdateTools();
}

void ToolWindowEdit::OnEditDrawRangeLow( wxCommandEvent& event )
{
  double value;
  if ( m_editDrawRangeLow->GetValue().ToDouble( &value ) )
  {
    BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
    double* range = bp->GetDrawRange();
    bp->SetDrawRange( value, range[1] );
    UpdateTools();
  }
}

void ToolWindowEdit::OnEditDrawRangeHigh( wxCommandEvent& event )
{
  double value;
  if ( m_editDrawRangeHigh->GetValue().ToDouble( &value ) )
  {
    BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
    double* range = bp->GetDrawRange();
    bp->SetDrawRange( range[0], value );
    UpdateTools();
  }
}

void ToolWindowEdit::OnEditExcludeRangeLow( wxCommandEvent& event )
{
  double value;
  if ( m_editExcludeRangeLow->GetValue().ToDouble( &value ) )
  {
    BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
    double* range = bp->GetExcludeRange();
    bp->SetExcludeRange( value, range[1] );
    UpdateTools();
  }
}

void ToolWindowEdit::OnEditExcludeRangeHigh( wxCommandEvent& event )
{
  double value;
  if ( m_editExcludeRangeHigh->GetValue().ToDouble( &value ) )
  {
    BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
    double* range = bp->GetExcludeRange();
    bp->SetExcludeRange( range[0], value );
    UpdateTools();
  }
}

void ToolWindowEdit::OnCheckSmooth( wxCommandEvent& event )
{
  for ( int i = 0; i < 3; i++ )
  {
    RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( i );
    Contour2D* c2d = view->GetContour2D();
    c2d->SetSmooth( event.IsChecked() );
  }
  UpdateTools();
}

void ToolWindowEdit::OnEditSmoothSD( wxCommandEvent& event )
{
  double value;

  if (m_editSmoothSD != NULL)
  {
    if ( m_editSmoothSD->GetValue().ToDouble( &value ) && value > 0 )
    {
      for ( int i = 0; i < 3; i++ )
      {
        RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( i );
        Contour2D* c2d = view->GetContour2D();
        c2d->SetSmoothSD( value );
      }
      UpdateTools();
    }
  }
}

void ToolWindowEdit::OnEditContourValue( wxCommandEvent& event )
{
  double value;
  if ( m_editContourValue->GetValue().ToDouble( &value ) && value > 0 )
  {
    BrushProperty* bp = MainWindow::GetMainWindowPointer()->GetBrushProperty();
    LayerMRI* mri = (LayerMRI*)bp->GetReferenceLayer();
    for ( int i = 0; i < 3; i++ )
    {
      RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( i );
      Contour2D* c2d = view->GetContour2D();
      if ( c2d->GetInputImage() )
        c2d->SetContourValue( value );
      else if ( mri )
      {
        c2d->SetInput( mri->GetSliceImageData( view->GetViewPlane() ), value, mri->GetSlicePosition()[i], mri->GetActiveFrame() ); 
        c2d->SetVisible( true );
      }
    }
    UpdateTools();
  }
}

void ToolWindowEdit::OnColorContour( wxColourPickerEvent& event )
{
  wxColour c = event.GetColour();
  for ( int i = 0; i < 3; i++ )
  {
    RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( i );
    Contour2D* c2d = view->GetContour2D();
    c2d->SetContourColor( c.Red()/255.0, c.Green()/255.0, c.Blue()/255.0 );
  }
}

