/**
 * @file  DialogCropVolume.h
 * @brief Dialog to crop volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.7 $
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

#include "DialogCropVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/config.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/file.h>
#include <wx/spinctrl.h> 
#include "MainWindow.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "VolumeCropper.h"
#include "vtkImageData.h"
#include "RenderView3D.h"

BEGIN_EVENT_TABLE( DialogCropVolume, wxFrame )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CLOSE" ),     DialogCropVolume::OnButtonClose )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_RESET" ),     DialogCropVolume::OnButtonReset )  
  EVT_BUTTON    ( XRCID( "ID_BUTTON_APPLY" ),     DialogCropVolume::OnButtonApply )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE_AS" ),   DialogCropVolume::OnButtonSaveAs ) 
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_X_MIN" ),       DialogCropVolume::OnSpinBound ) 
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_X_MAX" ),       DialogCropVolume::OnSpinBound ) 
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_Y_MIN" ),       DialogCropVolume::OnSpinBound ) 
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_Y_MAX" ),       DialogCropVolume::OnSpinBound ) 
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_Z_MIN" ),       DialogCropVolume::OnSpinBound )   
  EVT_SPINCTRL  ( XRCID( "ID_SPIN_Z_MAX" ),       DialogCropVolume::OnSpinBound )
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_X_MIN" ),       DialogCropVolume::OnSpinBoundText ) 
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_X_MAX" ),       DialogCropVolume::OnSpinBoundText ) 
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_Y_MIN" ),       DialogCropVolume::OnSpinBoundText ) 
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_Y_MAX" ),       DialogCropVolume::OnSpinBoundText ) 
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_Z_MIN" ),       DialogCropVolume::OnSpinBoundText ) 
  EVT_TEXT_ENTER      ( XRCID( "ID_SPIN_Z_MAX" ),       DialogCropVolume::OnSpinBoundText )
  EVT_SHOW     ( DialogCropVolume::OnShow )
  EVT_CLOSE    ( DialogCropVolume::OnClose )
END_EVENT_TABLE()


DialogCropVolume::DialogCropVolume( wxWindow* parent, LayerMRI* mri ) : Listener( "DialogCropVolume" )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_TOOLWINDOW_CROP_VOLUME") );
  m_spinRange[0]  = XRCCTRL( *this, "ID_SPIN_X_MIN", wxSpinCtrl );
  m_spinRange[1]  = XRCCTRL( *this, "ID_SPIN_X_MAX", wxSpinCtrl );
  m_spinRange[2]  = XRCCTRL( *this, "ID_SPIN_Y_MIN", wxSpinCtrl );
  m_spinRange[3]  = XRCCTRL( *this, "ID_SPIN_Y_MAX", wxSpinCtrl );
  m_spinRange[4]  = XRCCTRL( *this, "ID_SPIN_Z_MIN", wxSpinCtrl );
  m_spinRange[5]  = XRCCTRL( *this, "ID_SPIN_Z_MAX", wxSpinCtrl );
  
  SetVolume( mri );
}

DialogCropVolume::~DialogCropVolume()
{}

void DialogCropVolume::SetVolume( LayerMRI* mri )
{
  m_mri = mri;
  if ( m_mri )
  {
    int* dim = m_mri->GetImageData()->GetDimensions();
    for ( int i = 0; i < 6; i++ )
      m_spinRange[i]->SetRange( 0, dim[i/2]-1 );
  }
}

void DialogCropVolume::OnClose( wxCloseEvent& event )
{
  Hide();
}

void DialogCropVolume::OnShow( wxShowEvent& event )
{
  static bool bShowFrames = true;
  MainWindow* mainwnd = MainWindow::GetMainWindowPointer();
  RenderView3D* view = (RenderView3D*)mainwnd->GetRenderView( 3 );

//#if wxCHECK_VERSION(2,9,0)
#if wxVERSION_NUMBER > 2900
  if ( event.IsShown() )
#else
  if ( event.GetShow() )
#endif
  {
    if ( view )
    {
      bShowFrames = view->GetShowSliceFrames();
      view->SetShowSliceFrames( false );
    }
    
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x = config->Read( _T("/ToolWindowCropVolume/PosX"), 0L );
      int y = config->Read( _T("/ToolWindowCropVolume/PosY"), 0L );
      if ( x == 0 && y == 0 )
        Center();
      else
        Move( x, y );
    }
  }
  else
  {
    if ( view )
      view->SetShowSliceFrames( bShowFrames );
    if ( mainwnd->GetVolumeCropper() )
      mainwnd->GetVolumeCropper()->Show( false );
    if ( mainwnd->IsShown() )
      mainwnd->SetMode( 0 );
    
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x, y;
      GetPosition( &x, &y );
      config->Write( _T("/ToolWindowCropVolume/PosX"), (long) x );
      config->Write( _T("/ToolWindowCropVolume/PosY"), (long) y );
    }
  }
  
  event.Skip();
}

void DialogCropVolume::OnButtonClose( wxCommandEvent& event )
{
  Hide();
}

void DialogCropVolume::OnButtonReset( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->GetVolumeCropper()->Reset();
}

void DialogCropVolume::OnButtonSaveAs( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveVolumeAs();
}

void DialogCropVolume::OnButtonApply( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->GetVolumeCropper()->Apply();
  MainWindow::GetMainWindowPointer()->NeedRedraw();
}

void DialogCropVolume::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "LayerRemoved" && iData == m_mri )
  {
    VolumeCropper* vc = MainWindow::GetMainWindowPointer()->GetVolumeCropper();
    if ( vc )
      vc->SetEnabled( false );
    Hide();
  }
  else if ( iMessage == "CropBoundChanged" && iData == m_mri )
  {
    int* ext = MainWindow::GetMainWindowPointer()->GetVolumeCropper()->GetExtent();
    for ( int i = 0; i < 6; i++ )
      m_spinRange[i]->SetValue( ext[i] );
  }
}

void DialogCropVolume::OnSpinBound( wxSpinEvent& event )
{
  for ( int i = 0; i < 6; i++ )
  {
    if ( event.GetEventObject() == m_spinRange[i] )
    {
      MainWindow::GetMainWindowPointer()->GetVolumeCropper()->SetExtent( i, event.GetInt() );
      break;
    }
  }
}

void DialogCropVolume::OnSpinBoundText( wxCommandEvent& event )
{
  for ( int i = 0; i < 6; i++ )
  {
    if ( event.GetEventObject() == m_spinRange[i] )
    {
      long nValue = -1;
      event.GetString().ToLong( &nValue );
      if ( nValue <= m_spinRange[i]->GetMax() && nValue >=  m_spinRange[i]->GetMin() )
        MainWindow::GetMainWindowPointer()->GetVolumeCropper()->SetExtent( i, nValue );
      break;
    }
  }
}
