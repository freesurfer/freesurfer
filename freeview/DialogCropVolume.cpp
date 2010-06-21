/**
 * @file  DialogCropVolume.h
 * @brief Dialog to crop volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/06/21 18:37:50 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008-2009,
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

#include "DialogCropVolume.h"
#include <wx/xrc/xmlres.h>
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

BEGIN_EVENT_TABLE( DialogCropVolume, wxDialog )
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
END_EVENT_TABLE()


DialogCropVolume::DialogCropVolume( wxWindow* parent, LayerMRI* mri ) : Listener( "DialogCropVolume" )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_CROP_VOLUME") );
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

void DialogCropVolume::OnShow( wxShowEvent& event )
{
  static bool bShowFrames = true;
  MainWindow* mainwnd = MainWindow::GetMainWindowPointer();
  RenderView3D* view = (RenderView3D*)mainwnd->GetRenderView( 3 );
  if ( event.GetShow() )
  {
    bShowFrames = view->GetShowSliceFrames();
    view->SetShowSliceFrames( false );
  }
  else
  {
    view->SetShowSliceFrames( bShowFrames );
    if ( mainwnd->GetVolumeCropper() )
      mainwnd->GetVolumeCropper()->Show( false );
  }
  
  event.Skip();
}

void DialogCropVolume::OnButtonClose( wxCommandEvent& event )
{
  Close();
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
