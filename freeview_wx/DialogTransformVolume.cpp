/**
 * @file  DialogTransformVolume.h
 * @brief Dialog to transform volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
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



#include "DialogTransformVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/scrolbar.h>
#include <wx/notebook.h>
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "CommonDataStruct.h"

#if wxVERSION_NUMBER > 2900
#include <wx/bookctrl.h>
#endif

extern "C"
{
#include "mri.h"
}

BEGIN_EVENT_TABLE( DialogTransformVolume, wxDialog )
    EVT_BUTTON      ( XRCID( "ID_BUTTON_APPLY" ),     DialogTransformVolume::OnApply )
    EVT_BUTTON      ( XRCID( "ID_BUTTON_RESTORE" ),   DialogTransformVolume::OnRestore )
    EVT_BUTTON      ( XRCID( "ID_BUTTON_SAVEAS" ),    DialogTransformVolume::OnSaveAs )
    EVT_BUTTON      ( XRCID( "ID_BUTTON_SAVE_REG" ),  DialogTransformVolume::OnSaveReg )
    EVT_CHECKBOX    ( XRCID( "ID_CHECK_1" ),          DialogTransformVolume::OnCheck1 )
    EVT_CHECKBOX    ( XRCID( "ID_CHECK_2" ),          DialogTransformVolume::OnCheck2 )
    EVT_CHECKBOX    ( XRCID( "ID_CHECK_3" ),          DialogTransformVolume::OnCheck3 )
    EVT_TEXT        ( XRCID( "ID_TEXT_X" ),           DialogTransformVolume::OnTextTranslateX )
    EVT_TEXT        ( XRCID( "ID_TEXT_Y" ),           DialogTransformVolume::OnTextTranslateY )
    EVT_TEXT        ( XRCID( "ID_TEXT_Z" ),           DialogTransformVolume::OnTextTranslateZ )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_X" ),    DialogTransformVolume::OnScrollTranslateX )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_Y" ),    DialogTransformVolume::OnScrollTranslateY )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_Z" ),    DialogTransformVolume::OnScrollTranslateZ )    
    EVT_TEXT        ( XRCID( "ID_TEXT_SCALE_X" ),         DialogTransformVolume::OnTextScaleX )
    EVT_TEXT        ( XRCID( "ID_TEXT_SCALE_Y" ),         DialogTransformVolume::OnTextScaleY )
    EVT_TEXT        ( XRCID( "ID_TEXT_SCALE_Z" ),         DialogTransformVolume::OnTextScaleZ )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_SCALE_X" ),    DialogTransformVolume::OnScrollScaleX )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_SCALE_Y" ),    DialogTransformVolume::OnScrollScaleY )
    EVT_COMMAND_SCROLL        ( XRCID( "ID_SCROLLBAR_SCALE_Z" ),    DialogTransformVolume::OnScrollScaleZ )   
    EVT_NOTEBOOK_PAGE_CHANGED ( XRCID( "ID_NOTEBOOK" ),       DialogTransformVolume::OnPageChanged )
END_EVENT_TABLE()


DialogTransformVolume::DialogTransformVolume( wxWindow* parent ) : Listener( "DialogTransformVolume" )
{
  wxXmlResource::Get()->LoadDialog( this, 
				    parent, 
				    wxT("ID_DIALOG_TRANSFORM_VOLUME") );
  m_check[0] = XRCCTRL( *this, "ID_CHECK_1", wxCheckBox );
  m_check[1] = XRCCTRL( *this, "ID_CHECK_2", wxCheckBox );
  m_check[2] = XRCCTRL( *this, "ID_CHECK_3", wxCheckBox );
  m_check[0]->SetValue( true );
  m_choice[0] = XRCCTRL( *this, "ID_CHOICE_1", wxChoice );
  m_choice[1] = XRCCTRL( *this, "ID_CHOICE_2", wxChoice );
  m_choice[2] = XRCCTRL( *this, "ID_CHOICE_3", wxChoice );
  m_textAngle[0]  = XRCCTRL( *this, "ID_TEXT_ANGLE_1", wxTextCtrl );
  m_textAngle[1]  = XRCCTRL( *this, "ID_TEXT_ANGLE_2", wxTextCtrl );
  m_textAngle[2]  = XRCCTRL( *this, "ID_TEXT_ANGLE_3", wxTextCtrl );
  m_btnRestoreOriginal  = XRCCTRL( *this, "ID_BUTTON_RESTORE", wxButton );
  m_btnSaveAs     = XRCCTRL( *this, "ID_BUTTON_SAVEAS", wxButton );
  m_btnSaveReg    = XRCCTRL( *this, "ID_BUTTON_SAVE_REG", wxButton );
  
  m_radioAroundCenter = XRCCTRL( *this, "ID_RADIO_AROUND_CENTER", wxRadioButton );
  m_radioAroundCursor = XRCCTRL( *this, "ID_RADIO_AROUND_CURSOR", wxRadioButton );
  m_radioNearest      = XRCCTRL( *this, "ID_RADIO_NEAREST",       wxRadioButton );
  m_radioTrilinear    = XRCCTRL( *this, "ID_RADIO_TRILINEAR",     wxRadioButton );
  m_radioSinc         = XRCCTRL( *this, "ID_RADIO_SINC",          wxRadioButton );
  m_radioActiveVolume = XRCCTRL( *this, "ID_RADIO_ACTIVE_LAYER",  wxRadioButton );
  m_radioAllVolumes   = XRCCTRL( *this, "ID_RADIO_ALL_LAYERS",    wxRadioButton );
  
  m_notebook    = XRCCTRL( *this, "ID_NOTEBOOK", wxNotebook );
  m_textTranslate[0] =  XRCCTRL( *this, "ID_TEXT_X", wxTextCtrl );
  m_textTranslate[1] =  XRCCTRL( *this, "ID_TEXT_Y", wxTextCtrl );
  m_textTranslate[2] =  XRCCTRL( *this, "ID_TEXT_Z", wxTextCtrl );
  m_scrollTranslate[0] =  XRCCTRL( *this, "ID_SCROLLBAR_X", wxScrollBar );
  m_scrollTranslate[1] =  XRCCTRL( *this, "ID_SCROLLBAR_Y", wxScrollBar );
  m_scrollTranslate[2] =  XRCCTRL( *this, "ID_SCROLLBAR_Z", wxScrollBar );
  m_textScale[0] =  XRCCTRL( *this, "ID_TEXT_SCALE_X", wxTextCtrl );
  m_textScale[1] =  XRCCTRL( *this, "ID_TEXT_SCALE_Y", wxTextCtrl );
  m_textScale[2] =  XRCCTRL( *this, "ID_TEXT_SCALE_Z", wxTextCtrl );
  m_scrollScale[0] =  XRCCTRL( *this, "ID_SCROLLBAR_SCALE_X", wxScrollBar );
  m_scrollScale[1] =  XRCCTRL( *this, "ID_SCROLLBAR_SCALE_Y", wxScrollBar );
  m_scrollScale[2] =  XRCCTRL( *this, "ID_SCROLLBAR_SCALE_Z", wxScrollBar );
  
  UpdateUI();
}

DialogTransformVolume::~DialogTransformVolume()
{}

bool DialogTransformVolume::GetRotation( int nIndex_in, 
				      int& plane_out, 
				      double& angle_out )
{
  if ( nIndex_in < 0 || 
       nIndex_in > 2 || 
       !m_check[ nIndex_in ]->IsChecked() )
    return false;

  plane_out = m_choice[ nIndex_in ]->GetSelection();
  return m_textAngle[ nIndex_in ]->GetValue().ToDouble( &angle_out );
}

void DialogTransformVolume::OnApply( wxCommandEvent& event )
{
  int plane;
  double angle;
  if ( !m_check[0]->IsChecked() && 
       !m_check[1]->IsChecked() && 
       !m_check[2]->IsChecked() )
  {
    wxMessageDialog dlg( this, 
			 _("Must at least select one rotation."), 
			 _("Error"), 
			 wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  else if ( ( m_check[0]->IsChecked() && !GetRotation( 0, plane, angle ) ) ||
            ( m_check[1]->IsChecked() && !GetRotation( 1, plane, angle ) ) ||
            ( m_check[2]->IsChecked() && !GetRotation( 2, plane, angle ) ) )
  {
    wxMessageDialog dlg( this, 
			 _("Please enter correct rotation angle."), 
			 _("Error"), 
			 wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  DoRotate();
  
  UpdateUI();
}

void DialogTransformVolume::OnSaveAs( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveVolumeAs();
}

void DialogTransformVolume::OnSaveReg( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveRegistrationAs();
}

void DialogTransformVolume::OnRestore( wxCommandEvent& event )
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    layer->Restore();
    UpdateUI();
  }
}

void DialogTransformVolume::DoRotate()
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    RotationElement re;
    if ( m_radioAroundCursor->GetValue() )
    {
      MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->
        GetSlicePosition( re.Point );
      layer->RemapPositionToRealRAS( re.Point, re.Point );
    }
    else
    {
    // use center of the volume to rotate
      layer->GetRASCenter( re.Point );
    }
    re.SampleMethod = SAMPLE_TRILINEAR;
    if ( m_radioNearest->GetValue() )
      re.SampleMethod = SAMPLE_NEAREST;
    else if ( m_radioSinc->GetValue() )
      re.SampleMethod = SAMPLE_SINC;
    
    std::vector<RotationElement> rotations;
    for ( int i = 0; i < 3; i++ )
    {
      if ( GetRotation( i, re.Plane, re.Angle ) )
      {
        rotations.push_back( re );
      }
    }
    MainWindow::GetMainWindowPointer()->RotateVolume( rotations, m_radioAllVolumes->GetValue() );
  }
}

void DialogTransformVolume::OnCheck1( wxCommandEvent& event )
{
  m_choice[0]->Enable( event.IsChecked() );
  m_textAngle[0]->Enable( event.IsChecked() );
}


void DialogTransformVolume::OnCheck2( wxCommandEvent& event )
{
  m_choice[1]->Enable( event.IsChecked() );
  m_textAngle[1]->Enable( event.IsChecked() );
}


void DialogTransformVolume::OnCheck3( wxCommandEvent& event )
{
  m_choice[2]->Enable( event.IsChecked() );
  m_textAngle[2]->Enable( event.IsChecked() );
}

void DialogTransformVolume::OnTextTranslateX( wxCommandEvent& event )
{
  RespondTextTranslate( 0 ); 
}

void DialogTransformVolume::OnTextTranslateY( wxCommandEvent& event )
{
  RespondTextTranslate( 1 ); 
}

void DialogTransformVolume::OnTextTranslateZ( wxCommandEvent& event )
{
  RespondTextTranslate( 2 ); 
}

void DialogTransformVolume::OnScrollTranslateX( wxScrollEvent& event )
{
  RespondScrollTranslate( 0 );
}

void DialogTransformVolume::OnScrollTranslateY( wxScrollEvent& event )
{
  RespondScrollTranslate( 1 );
}

void DialogTransformVolume::OnScrollTranslateZ( wxScrollEvent& event )
{
  RespondScrollTranslate( 2 );
}

void DialogTransformVolume::RespondTextTranslate( int n )
{
  if ( IsShown() )
  {
    double dvalue;
    if ( m_textTranslate[n]->GetValue().ToDouble( &dvalue ) )
    {
      LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
      if ( layer )
      {
        double pos[3];
        layer->GetTranslate( pos );
        pos[n] = dvalue;
        layer->Translate( pos );
        MainWindow::GetMainWindowPointer()->NeedRedraw();
        
        double* vs = layer->GetWorldVoxelSize();
        int range = m_scrollTranslate[n]->GetRange();
        m_scrollTranslate[n]->SetThumbPosition( range/2 + (int)( pos[n] / vs[n] ) );
        m_btnRestoreOriginal->Enable();
        UpdateUI( 1 );
      }
    }
  }
}

void DialogTransformVolume::RespondScrollTranslate( int n )
{
  if ( IsShown() )
  {
    LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
    if ( layer )
    {
      double pos[3];
      layer->GetTranslate( pos );
      int range = m_scrollTranslate[n]->GetRange();
      int npos = m_scrollTranslate[n]->GetThumbPosition();
      double* vs = layer->GetWorldVoxelSize();
      pos[n] = ( npos - range/2 ) * vs[n];    
      layer->Translate( pos );
      MainWindow::GetMainWindowPointer()->NeedRedraw();
      
      m_textTranslate[n]->ChangeValue( wxString() << pos[n] );
      m_btnRestoreOriginal->Enable();
      UpdateUI( 1 );
    }
  }
}


void DialogTransformVolume::OnTextScaleX( wxCommandEvent& event )
{
  RespondTextScale( 0 ); 
}

void DialogTransformVolume::OnTextScaleY( wxCommandEvent& event )
{
  RespondTextScale( 1 ); 
}

void DialogTransformVolume::OnTextScaleZ( wxCommandEvent& event )
{
  RespondTextScale( 2 ); 
}

void DialogTransformVolume::OnScrollScaleX( wxScrollEvent& event )
{
  RespondScrollScale( 0 );
}

void DialogTransformVolume::OnScrollScaleY( wxScrollEvent& event )
{
  RespondScrollScale( 1 );
}

void DialogTransformVolume::OnScrollScaleZ( wxScrollEvent& event )
{
  RespondScrollScale( 2 );
}

void DialogTransformVolume::RespondTextScale( int n )
{
  if ( IsShown() )
  {
    double dvalue;
    if ( m_textScale[n]->GetValue().ToDouble( &dvalue ) && dvalue > 0 )
    {
      LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
      if ( layer )
      {
        double scale[3];
        layer->GetScale( scale );
        scale[n] = dvalue;
        layer->Scale( scale );
        MainWindow::GetMainWindowPointer()->NeedRedraw();
        
        if ( dvalue >= 1 )
          m_scrollScale[n]->SetThumbPosition( 50 + (int)( (dvalue-1.0)*50 ) );
        else        
          m_scrollScale[n]->SetThumbPosition( 50 - (int)( (1.0-dvalue)*100 ) );
          
        m_btnRestoreOriginal->Enable();
        UpdateUI( 0 );
      }
    }
  }
}

void DialogTransformVolume::RespondScrollScale( int n )
{
  LayerMRI* layer = ( LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    double scale[3];
    layer->GetScale( scale );
    int npos = m_scrollScale[n]->GetThumbPosition();
    if ( npos >= 50 )
      scale[n] = ( npos - 50 ) / 50.0 + 1.0;
    else 
      scale[n] = ( npos - 50 ) / 100.0 + 1.0;    
    layer->Scale( scale );
    MainWindow::GetMainWindowPointer()->NeedRedraw();
    
    m_textScale[n]->ChangeValue( wxString() << scale[n] );
    m_btnRestoreOriginal->Enable();
    UpdateUI( 0 );
  }
}

//#if wxCHECK_VERSION( 2, 9, 0 )
#if wxVERSION_NUMBER > 2900
void DialogTransformVolume::OnPageChanged( wxBookCtrlEvent& event )
#else
void DialogTransformVolume::OnPageChanged( wxNotebookEvent& event )
#endif
{
}

void DialogTransformVolume::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ActiveLayerChanged" )
  {
    if ( IsVisible() )
      UpdateUI();
  }
}

// scope: 0 => translate related, 1 => scale related, 2 => both
void DialogTransformVolume::UpdateUI( int scope )
{
  LayerMRI* layer = (LayerMRI* )MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  if ( layer )
  {
    if ( scope == 0 || scope == 2 )
    {
      double* vs = layer->GetWorldVoxelSize();
      double* ws = layer->GetWorldSize();
      double pos[3];
      layer->GetTranslate( pos );
      for ( int i = 0; i < 3; i++ )
      {
        int range = (int)( ws[i] / vs[i] + 0.5 ) * 2;
        int npos = (int)(pos[i] / vs[i]) + range/2;
        m_scrollTranslate[i]->SetScrollbar( npos, 0, range, 2, true );
        m_textTranslate[i]->ChangeValue( wxString() << pos[i] );
      }
    }
    if ( scope == 1 || scope == 2 )
    {
      double scale[3];
      layer->GetScale( scale );
      for ( int i = 0; i < 3; i++ )
      {
        if ( scale[i] >= 1 )
          m_scrollScale[i]->SetThumbPosition( 50 + (int)( (scale[i]-1.0)*50 ) );
        else        
          m_scrollScale[i]->SetThumbPosition( 50 - (int)( (1.0-scale[i])*100 ) );
          
        m_textScale[i]->ChangeValue( wxString() << scale[i] );
      }
    }
    
    m_btnRestoreOriginal->Enable( layer->IsTransformed() );
    m_btnSaveReg->Enable( layer->IsTransformed() );
    m_btnSaveAs->Enable( layer->IsTransformed() );
  }
}
