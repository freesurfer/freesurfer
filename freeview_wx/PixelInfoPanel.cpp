/**
 * @file  PixelInfoPanel.h
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

#include "PixelInfoPanel.h"
#include <wx/xrc/xmlres.h>
#include <wx/aui/auibook.h>
#include <wx/listctrl.h>
#include <wx/textctrl.h>
#include <wx/config.h>
#include "stdlib.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPLabel.h"
#include "LayerSurface.h"
#include "SurfaceOverlay.h"
#include "SurfaceAnnotation.h"
#include "MyUtils.h"
#include "LayerProperties.h"

#define  DEFAULT_NAME_LENGTH    12

BEGIN_EVENT_TABLE(PixelInfoListCtrl, wxListCtrl)
//    EVT_LIST_BEGIN_DRAG(PIXEL_INFO_LISTCTRL_ID, PixelInfoListCtrl::OnBeginDrag)
  EVT_LIST_COL_BEGIN_DRAG   ( PIXEL_INFO_LISTCTRL_ID,         PixelInfoListCtrl::OnDrag )
  EVT_LIST_COL_DRAGGING     ( PIXEL_INFO_LISTCTRL_ID,         PixelInfoListCtrl::OnDrag )
  EVT_LIST_COL_END_DRAG     ( PIXEL_INFO_LISTCTRL_ID,         PixelInfoListCtrl::OnDrag )
  EVT_TEXT_ENTER            ( PIXEL_INFO_INLINE_EDITOR_ID,    PixelInfoListCtrl::OnFinishEditing )
  EVT_MENU                  ( MENU_SHOW_VOXEL_COORDINATES,    PixelInfoListCtrl::OnShowVoxelCoordinates )
  EVT_MENU                  ( MENU_SHOW_SURFACE_COORDINATES,  PixelInfoListCtrl::OnShowSurfaceCoordinates )
  EVT_MENU                  ( MENU_SHOW_SURFACE_CURVATURE,    PixelInfoListCtrl::OnShowSurfaceCurvature )
  EVT_MENU                  ( MENU_SHOW_SHORT_NAME,           PixelInfoListCtrl::OnShowShortName )
  EVT_SIZE                  ( PixelInfoListCtrl::OnSize )
  EVT_LEFT_DOWN             ( PixelInfoListCtrl::OnButtonDown )
  EVT_RIGHT_DOWN            ( PixelInfoListCtrl::OnButtonDown )
  EVT_CONTEXT_MENU          ( PixelInfoListCtrl::OnContextMenu )
END_EVENT_TABLE()

PixelInfoListCtrl::PixelInfoListCtrl( wxWindow* parent, const wxString& name ) :
    wxListCtrl( parent,
                PIXEL_INFO_LISTCTRL_ID,
                wxDefaultPosition,
                wxDefaultSize,
                (wxLC_REPORT | wxLC_SINGLE_SEL | wxLC_VIRTUAL | wxBORDER_SUNKEN ) &~wxHSCROLL ),
    Listener( "PixelInfoCtrl" ),
    m_attr(*wxBLACK, wxColour(230, 230, 230), wxNullFont)
{
  InsertColumn( 0, name );
  InsertColumn( 1, _("") );

  for ( int i = 0; i < 3; i++ )
    m_dRASPosition[i] = 0;

  m_textEditor = new wxTextCtrl( this, PIXEL_INFO_INLINE_EDITOR_ID, _(""),
                                 wxDefaultPosition,
                                 wxDefaultSize,
                                 wxTE_PROCESS_ENTER );
  m_textEditor->Hide();
  m_bShowVoxelCoordinates = true;
  m_bShowSurfaceCoordinates = true;
  m_bShowSurfaceCurvature = true;
  m_bShowShortName = false;

  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int w = config->Read( _T("/VoxelInfo/ColumnWidth/") + name, 100L );
    SetColumnWidth( 0, w );
  }
  else
    SetColumnWidth( 0, 100 );
  
  Reset();
}

PixelInfoListCtrl::~PixelInfoListCtrl()
{
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    wxListItem item;
    GetColumn( 0, item );
    wxString name = item.GetText();
    config->Write( _T("/VoxelInfo/ColumnWidth/") + name, (long)GetColumnWidth( 0 ) );
  } 
}

wxString PixelInfoListCtrl::OnGetItemText(long item, long column) const
{
  wxString ret;
  if ( item < (int)m_listName.Count() )
  {
    if ( column == 0 )
      ret = m_listName[item];
    else if ( column == 1 )
      ret = m_listValue[item];
  }
  return ret;
}

wxListItemAttr* PixelInfoListCtrl::OnGetItemAttr(long item) const
{
  return item % 2 ? NULL : (wxListItemAttr *)&m_attr;
}

void PixelInfoListCtrl::OnSize( wxSizeEvent& event )
{
  event.Skip();
  ResizeColumns();
}

void PixelInfoListCtrl::OnDrag( wxListEvent& event )
{
  event.Skip();
  ResizeColumns();
}

void PixelInfoListCtrl::OnButtonDown( wxMouseEvent& event )
{
  event.Skip();

  wxPoint pt = event.GetPosition();
  int nStyle = wxLIST_HITTEST_ONITEM;
  long nCol;
  long nRow = HitTest( pt, nStyle, &nCol );
  if ( event.LeftDown() && ( m_bShowVoxelCoordinates && nRow >= 0 && nRow < GetItemCount() ) )
  {
    if ( m_listPtr[nRow] != 0 )
    {
      wxRect rect;
      GetItemRect( nRow, rect );
      rect.SetX( rect.GetX() + GetColumnWidth( 0 ) );
  //  std::cout << pt.x << " " << pt.y << "  " << rect.GetX() << " " << rect.GetY() << std::endl;
      rect.Inflate( 1, 1 );
      rect.SetWidth( rect.GetWidth() - GetColumnWidth( 0 ) - 18 );
  
      m_textEditor->Show();
      m_textEditor->SetFocus();
      m_textEditor->SetSize( rect );
  
      wxString strg = m_listValue[nRow];
      if ( strg.Find( _("Vertex") ) != wxNOT_FOUND )
      {
        strg = strg.Mid( 6 ).Trim( false );
        if ( strg.Find( _("[") ) != wxNOT_FOUND )
          strg = strg.Left( strg.Find( _("[") ) ).Trim( false ).Trim( true );
      }
      else if ( strg.Find( _("[") ) != wxNOT_FOUND && strg.Find( _("]") ) != wxNOT_FOUND )
      {
        strg = strg.Mid( strg.Find( _("[") ) + 1 );
        strg = strg.Left( strg.Find( _("]") ) );
      }
      m_textEditor->SetValue( strg );
      m_nRowEdited = nRow;
    }
  }
  else
  {
    m_textEditor->Hide();
  }

  if ( event.RightDown() )
  {
    wxMenu menu;
    wxMenuItem* item = menu.AppendCheckItem( MENU_SHOW_VOXEL_COORDINATES, _("Show Voxel Coordinate") );
    item->Check( m_bShowVoxelCoordinates );
    if ( !MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->IsEmpty() )
    {
      item = menu.AppendCheckItem( MENU_SHOW_SURFACE_COORDINATES, _("Show Surface Coordinate") );
      item->Check( m_bShowSurfaceCoordinates );
      item = menu.AppendCheckItem( MENU_SHOW_SURFACE_CURVATURE, _("Show Curvature Value") );
      item->Check( m_bShowSurfaceCurvature );
    }
    menu.AppendSeparator();
    item = menu.AppendCheckItem( MENU_SHOW_SHORT_NAME, _("Show Short Label Names") );
    item->Check( m_bShowShortName );
    PopupMenu( &menu, pt.x, pt.y );
  }
}


void PixelInfoListCtrl::ResizeColumns()
{
  wxSize sz = GetClientSize();
  wxRect rc = GetViewRect();
  SetColumnWidth( 1, sz.GetWidth() - GetColumnWidth( 0 ) - ( rc.GetHeight() > sz.GetHeight() ? 18 : 0 ) );
}

void PixelInfoListCtrl::AddItem( const wxString& name, 
                                 const wxString& value,
                                 long userData, 
                                 bool bUpdateCount )
{
  m_listName.Add( name );
  m_listValue.Add( value );
  m_listPtr.Add( userData );
  if ( bUpdateCount )
    SetItemCount( m_listName.Count() );
}

void PixelInfoListCtrl::SetRASPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
    m_dRASPosition[i] = pos[i];

  UpdateList();
}

void PixelInfoListCtrl::Reset()
{
  m_listName.Clear();
  m_listValue.Clear();
  m_listPtr.Clear();

//  AddItem( "RAS" );
}

void PixelInfoListCtrl::OnInternalIdle()
{
  if ( m_bNeedUpdate )
  {
    DoUpdateList();
    m_bNeedUpdate = false;
  }

  wxListCtrl::OnInternalIdle();
}

void PixelInfoListCtrl::UpdateList()
{
  m_bNeedUpdate = true;
}

void PixelInfoListCtrl::DoUpdateList()
{
  Reset();

  LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerCollection* lc_surf = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );
  if ( !lc_mri || !lc_surf )
  {
    // collection does not exist. program must be in exiting process.
    return;
  }
  if ( lc_mri->IsEmpty() && lc_surf->IsEmpty() )
  {
    SetItemCount( 0 );
    return;
  }
  else
    AddItem( _("RAS"), _(""), 1, false );
  
  if ( !lc_mri->IsEmpty() )
  {
    LayerMRI* mri = ( LayerMRI* )lc_mri->GetLayer( 0 );
    double ras[3];
    mri->RemapPositionToRealRAS( m_dRASPosition, ras );
    m_listValue[0].Printf( _("%.2f, %.2f, %.2f"), ras[0], ras[1], ras[2] );

    int nIndex[3];
    std::vector<Layer*> layers = lc_mri->GetLayers();
    for ( size_t i = 0; i < layers.size(); i++ )
    {      
      if ( layers[i]->GetProperties()->GetShowInfo() )
      {
        ((LayerMRI*)layers[i])->RASToOriginalIndex( ras, nIndex );
  //      double dvalue = ( (LayerMRI*)layers[i] )->GetVoxelValueByOriginalIndex( nIndex[0], nIndex[1], nIndex[2] );
        double dvalue = ( (LayerMRI*)layers[i] )->GetVoxelValue( m_dRASPosition );  
        wxString coordStrg;
        if ( m_bShowVoxelCoordinates )
          coordStrg.Printf( _("[%d, %d, %d] "), nIndex[0], nIndex[1], nIndex[2] );
        
        wxString labelStrg;
        if (layers[i]->IsTypeOf("PLabel"))
        {
          labelStrg = wxString::FromAscii( ( (LayerPLabel*)layers[i] )->GetLabelName( m_dRASPosition ).c_str() );
        }
        else
        {
          labelStrg = wxString::FromAscii( ( (LayerMRI*)layers[i] )->GetLabelName( dvalue ).c_str() );
          if ( m_bShowShortName )
            labelStrg = GetShortName( labelStrg );
        }
        
        AddItem( wxString::FromAscii( layers[i]->GetName() ),
                AppendSpaceString( ( wxString() << dvalue ) ) + coordStrg + _("  ") + labelStrg,
                (long)layers[i],
                false );
      }
    }
  }
  if ( !lc_surf->IsEmpty() )
  {
    if ( lc_mri->IsEmpty() )
    {
      double ras[3] = { m_dRASPosition[0], m_dRASPosition[1], m_dRASPosition[2] };
      // mri->RemapPositionToRealRAS( m_dRASPosition, ras );
      m_listValue[0].Printf( _("%.2f, %.2f, %.2f"), ras[0], ras[1], ras[2] );
    }

    std::vector<Layer*> layers = lc_surf->GetLayers();
    for ( size_t i = 0; i < layers.size(); i++ )
    {
      if ( layers[i]->GetProperties()->GetShowInfo() )
      {
        LayerSurface* surf = ( LayerSurface* )layers[i];
        AddSurfaceItem( surf, m_dRASPosition, false );
      }
    }
  }

  SetItemCount( m_listName.Count() );
}

void PixelInfoListCtrl::AddSurfaceItem( LayerSurface* surf, double* pos, bool bUpdateCount )
{
  if ( m_bShowSurfaceCoordinates )
  {
    double sf_pos[3];
    surf->GetSurfaceRASAtTarget( pos, sf_pos );
    wxString strg;
    strg.Printf( _("[%.2f, %.2f, %.2f]"), sf_pos[0], sf_pos[1], sf_pos[2] );
    wxString text = ( AppendSpaceString( _("Coordinate") ) << strg );
    AddItem( wxString::FromAscii( surf->GetName() ), text, (long)surf );
  }
  int nVertex = surf->GetVertexIndexAtTarget( pos, NULL );
  if ( nVertex < 0 )
  {
    if ( !m_bShowSurfaceCoordinates )
      AddItem( wxString::FromAscii( surf->GetName() ), _(""), (long)surf, bUpdateCount );
  }
  else
  {
    wxString text = AppendSpaceString( _("Vertex") );
    text << nVertex;
    double ras[3];
    surf->GetSurfaceRASAtVertex( nVertex, ras ); 
    wxString strg;
    strg.Printf( _("[%.2f, %.2f, %.2f]"), ras[0], ras[1], ras[2] );
    text = ( AppendSpaceString(text) << strg );    
    AddItem( ( m_bShowSurfaceCoordinates ? wxString("") : wxString::FromAscii( surf->GetName() ) ), text, (long)surf );
    
    double vec[3];
    surf->GetNormalAtVertex( nVertex, vec ); 
    strg.Printf( _("[%.2f, %.2f, %.2f]"), vec[0], vec[1], vec[2] );
    text = ( AppendSpaceString( _("Normal") ) << strg );
    AddItem( _(""), text );
    
    if ( surf->GetActiveVector() >= 0 )
    {
      surf->GetVectorAtVertex( nVertex, vec ); 
      wxString strg;
      strg.Printf( _("[%.2f, %.2f, %.2f]"), vec[0], vec[1], vec[2] );
      text = ( AppendSpaceString( _("Vector") ) << strg );
      AddItem( _(""), text );
    }
    
    if ( surf->HasCurvature() && m_bShowSurfaceCurvature )
      AddItem( _(""), ( AppendSpaceString( _("Curvature") ) << surf->GetCurvatureValue( nVertex ) ), 0, bUpdateCount );
    
    
    int nOverlays = surf->GetNumberOfOverlays();
    for ( int i = 0; i < nOverlays; i++ )
    {
      SurfaceOverlay* overlay = surf->GetOverlay( i );
      wxString strg = overlay->GetName();
      AddItem( _(""), ( AppendSpaceString( strg ) << overlay->GetDataAtVertex( nVertex ) ), 0, bUpdateCount ); 
    }
    
    int nAnnotations = surf->GetNumberOfAnnotations();
    for ( int i = 0; i < nAnnotations; i++ )
    {
      SurfaceAnnotation* annot = surf->GetAnnotation( i );
      wxString strg = annot->GetName();
      AddItem( _(""), ( AppendSpaceString( strg ) << annot->GetAnnotationNameAtVertex( nVertex ).c_str() ), 0, bUpdateCount ); 
    }
  }
}

wxString PixelInfoListCtrl::GetShortName( wxString& name )
{
  if ( name.IsEmpty() || name.Find( _("-") ) == wxNOT_FOUND )
    return name;

  wxString short_name;
  wxArrayString subs = MyUtils::SplitString( name, _("-") );
  for ( int i = 0; i < (int)subs.Count(); i++ )
  {
    short_name += subs[i].Left( 1 );
  }
  return short_name;
}

void PixelInfoListCtrl::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || iMsg == "LayerMoved" ||
       iMsg == "MouseRASPositionChanged" || iMsg == "SurfaceOverlayAdded" ||
       iMsg == "SurfaceAnnotationAdded" )
  {
    HideEditor();
    UpdateList();
  }
}

void PixelInfoListCtrl::OnFinishEditing( wxCommandEvent& event )
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );

  // see if entered text is valid
  wxArrayString sa = MyUtils::SplitString( m_textEditor->GetValue(), _(",") );
  long ptr = m_listPtr[m_nRowEdited];
  Layer* layer_ptr = (Layer*)ptr;
  if ( ptr == 1 || ( layer_ptr && layer_ptr->IsTypeOf( "MRI" ) ) )
  {
    if ( sa.Count() < 3 )
      sa = MyUtils::SplitString( m_textEditor->GetValue(), _(" ") );

    if ( sa.Count() < 3 )
    {
      cerr << "Invalid coordinate string. Make sure they are three numbers." << endl;
      return;
    }
  }

  if ( ptr == 1 ) // RAS
  {
    double ras[3];
    if ( sa[0].ToDouble( ras ) && sa[1].ToDouble( ras+1 ) && sa[2].ToDouble( ras+2 ) )
    {
      wxListItem item;
      GetColumn( 0, item );
      LayerMRI* layer = (LayerMRI*)lc->GetLayer( 0 );
      if ( layer )
        layer->RASToTarget( ras, ras );
      if ( item.GetText() == _("Cursor") )
      {
        lc->SetCursorRASPosition( ras );
        lcm->SetSlicePosition( ras );
      }
      else if ( item.GetText() == _("Mouse") )
        lc->SetCurrentRASPosition( ras );
      UpdateList();
      m_textEditor->Hide();
    }
    else
    {
      cerr << "Invalid coordinate string. Make sure they are three numbers." << endl;
    }
  }
  else if ( layer_ptr && layer_ptr->IsTypeOf( "MRI" ) ) // voxel
  {
    long x, y, z;
    if ( sa.Count() < 3 )
    {
      cerr << "Invalid voxel coordinate string. Make sure they are three numbers." << endl;
      return;
    }
    int n = sa[0].Find( wxChar('['), true );
    if ( n != wxNOT_FOUND )
      sa[0] = sa[0].Mid( n+1 );
    n = sa[2].Find( wxChar(']') );
    if ( n != wxNOT_FOUND )
      sa[2] = sa[2].Left( n );
    if ( sa[0].ToLong( &x ) && sa[1].ToLong( &y ) && sa[2].ToLong( &z ) )
    {
      int nv[3] = { x, y, z };
      double ras[3];
      wxListItem item;
      GetColumn( 0, item );
      LayerMRI* layer = (LayerMRI*)layer_ptr;
      layer->OriginalIndexToRAS( nv, ras );
      layer->RASToTarget( ras, ras );
      if ( item.GetText() == _("Cursor") )
      {
        lc->SetCursorRASPosition( ras );
        lcm->SetSlicePosition( ras );
      }
      else if ( item.GetText() == _("Mouse") )
        lc->SetCurrentRASPosition( ras );
      UpdateList();
      m_textEditor->Hide();
    }
    else
    {
      cerr << "Invalid voxel coordinate string. Make sure they are three numbers." << endl;
    }
  }  
  else if ( layer_ptr && layer_ptr->IsTypeOf( "Surface" )  ) // surface
  {
    wxString strg = m_textEditor->GetValue();
    LayerSurface* layer = (LayerSurface*)layer_ptr;
    double ras[3];
    bool bSuccess = false;
    if ( m_listValue[m_nRowEdited].Find( _("Coord") ) == 0 )  // coordinate item
    {
      sa = MyUtils::SplitString( strg, _(",") );
      if ( sa.Count() < 3 )
        sa = MyUtils::SplitString( m_textEditor->GetValue(), _(" ") );
      if ( sa.Count() >= 3 && sa[0].ToDouble( ras ) && sa[1].ToDouble( ras+1 ) && sa[2].ToDouble( ras+2 ) )
      {
        layer->GetTargetAtSurfaceRAS( ras, ras );
        bSuccess = true;
      }
    }
    else        // vertex index item
    {
      long nIndex;
      if ( strg.ToLong( &nIndex ) && layer->GetTargetAtVertex( nIndex, ras ) )
        bSuccess = true;
    }
    if ( bSuccess )
    {
      wxListItem item;
      GetColumn( 0, item );
      if ( item.GetText() == _("Cursor") )
      {
        lc->SetCursorRASPosition( ras );
        lcm->SetSlicePosition( ras );
      }
      else if ( item.GetText() == _("Mouse") )
        lc->SetCurrentRASPosition( ras );
          
      UpdateList();
      m_textEditor->Hide();
    }
    else
      cerr << "Invalid index or coordinate string." << endl;
  }
}

void PixelInfoListCtrl::HideEditor()
{
  m_textEditor->Hide();
}

void PixelInfoListCtrl::SetShowVoxelCoordinates( bool bShow )
{
  m_bShowVoxelCoordinates = bShow;
  HideEditor();
  UpdateList();
}

void PixelInfoListCtrl::SetShowSurfaceCoordinates( bool bShow )
{
  m_bShowSurfaceCoordinates = bShow;
  HideEditor();
  UpdateList();
}

void PixelInfoListCtrl::SetShowSurfaceCurvature( bool bShow )
{
  m_bShowSurfaceCurvature = bShow;
  HideEditor();
  UpdateList();
}

void PixelInfoListCtrl::SetShowShortName( bool bShow )
{
  m_bShowShortName = bShow;
  HideEditor();
  UpdateList();
}

void PixelInfoListCtrl::OnContextMenu( wxContextMenuEvent& event )
{}

void PixelInfoListCtrl::OnShowVoxelCoordinates( wxCommandEvent& event )
{
  SetShowVoxelCoordinates( event.IsChecked() );
}

void PixelInfoListCtrl::OnShowSurfaceCoordinates( wxCommandEvent& event )
{
  SetShowSurfaceCoordinates( event.IsChecked() );
}

void PixelInfoListCtrl::OnShowSurfaceCurvature( wxCommandEvent& event )
{
  SetShowSurfaceCurvature( event.IsChecked() );
}

void PixelInfoListCtrl::OnShowShortName( wxCommandEvent& event )
{
  SetShowShortName( event.IsChecked() );
}

wxString PixelInfoListCtrl::AppendSpaceString( const wxString& input, int nDefaultLength )
{
  int nLen = nDefaultLength;
  if ( nLen <= 0 )
    nLen = DEFAULT_NAME_LENGTH;
  
  wxString str = input;
  for ( int i = 0; i < nLen - (int)input.Len(); i++ )
    str += _(" ");
  
  str += _(" \t");
  
  return str;
}

BEGIN_EVENT_TABLE( PixelInfoPanel, wxPanel )
/* EVT_IDLE(PixelInfoPanel::OnIdle)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_X")),   PixelInfoPanel::OnCheckBoxSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_Y")),   PixelInfoPanel::OnCheckBoxSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_Z")),   PixelInfoPanel::OnCheckBoxSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_X")),   PixelInfoPanel::OnEditSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_Y")),   PixelInfoPanel::OnEditSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_Z")),   PixelInfoPanel::OnEditSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_X")),   PixelInfoPanel::OnSliderSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_Y")),   PixelInfoPanel::OnSliderSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_Z")),   PixelInfoPanel::OnSliderSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_ATTACH_SLICE")),  PixelInfoPanel::OnCheckBoxAttachSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_TRANSPARENCY")), PixelInfoPanel::OnSliderTransparency)*/
END_EVENT_TABLE()

PixelInfoPanel::PixelInfoPanel( wxWindow* parent ) :
    wxPanel( parent, wxID_ANY),
    Broadcaster( "PixelInfoPanel"),
    Listener( "PixelInfoPanel" )
{
  wxBoxSizer* sizer = new wxBoxSizer( wxHORIZONTAL );

  m_listCtrlCursor = new PixelInfoListCtrl( this, _("Cursor") );
  m_listCtrlMouse = new PixelInfoListCtrl( this, _("Mouse") );

  sizer->Add( m_listCtrlCursor, 1, wxEXPAND );
  sizer->AddSpacer( 3 );
  sizer->Add( m_listCtrlMouse, 1, wxEXPAND );

  this->SetSizer( sizer );
  this->Layout();

  AddListener( m_listCtrlCursor );
  AddListener( m_listCtrlMouse );
}

PixelInfoPanel::~PixelInfoPanel()
{}

void PixelInfoPanel::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "MouseRASPositionChanged" )
  {
    m_listCtrlCursor->BlockListen( true );
    m_listCtrlMouse->SetRASPosition( MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetCurrentRASPosition() );
    m_listCtrlMouse->HideEditor();
  }
  else if ( iMsg == "CursorRASPositionChanged" )
  {
    m_listCtrlMouse->BlockListen( true );
    m_listCtrlCursor->SetRASPosition( MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetCursorRASPosition() );
    m_listCtrlCursor->HideEditor();
  }
  else if ( iMsg == "LayerActiveFrameChanged" || iMsg == "LayerShowInfoChanged" )
  {
    m_listCtrlMouse->UpdateList();
    m_listCtrlCursor->UpdateList();
  }

  SendBroadcast( iMsg, iData, sender );
  m_listCtrlMouse->BlockListen( false );
  m_listCtrlCursor->BlockListen( false );
}

void PixelInfoPanel::ToggleShowVoxelCoordinates()
{
  bool bShow = m_listCtrlCursor->GetShowVoxelCoordinates();
  m_listCtrlCursor->SetShowVoxelCoordinates( !bShow );
  m_listCtrlMouse->SetShowVoxelCoordinates( !bShow );
}
