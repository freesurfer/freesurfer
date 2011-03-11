/**
 * @file  DialogSavePoint.h
 * @brief Dialog to save volume.
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



#include "DialogSavePoint.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/file.h>
#include <wx/combobox.h>
#include "MainWindow.h"
#include "MyUtils.h"

BEGIN_EVENT_TABLE( DialogSavePoint, wxDialog )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CLOSE" ),     DialogSavePoint::OnClose )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE" ),      DialogSavePoint::OnSave )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_DELETE" ),    DialogSavePoint::OnDelete )   
  EVT_COMBOBOX  ( XRCID( "ID_COMBO_GROUP_NAME" ), DialogSavePoint::OnComboGroupName )  
  EVT_SHOW      ( DialogSavePoint::OnShow ) 
END_EVENT_TABLE()


DialogSavePoint::DialogSavePoint( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_SAVE_POINT") );
  m_comboGroupName  = XRCCTRL( *this, "ID_COMBO_GROUP_NAME", wxComboBox );
  m_comboPointName  = XRCCTRL( *this, "ID_COMBO_POINT_NAME", wxComboBox );
}

DialogSavePoint::~DialogSavePoint()
{}

wxString DialogSavePoint::GetGroupName()
{
  return m_comboGroupName->GetValue().Trim( true ).Trim( false );
}

wxString DialogSavePoint::GetPointName()
{
  return m_comboPointName->GetValue().Trim( true ).Trim( false );
}

void DialogSavePoint::OnClose( wxCommandEvent& event )
{
  Close();
}

void DialogSavePoint::OnSave( wxCommandEvent& event )
{
  if ( !ValidateInput() )
    return;
  
  wxString group_name = m_comboGroupName->GetValue().Trim( true ).Trim( false );
  wxString point_name = m_comboPointName->GetValue().Trim( true ).Trim( false );
  double ras[3];
  MainWindow::GetMainWindowPointer()->GetCursorRAS( ras );
  bool bFound = false;
  for ( size_t i = 0; i < m_points.size(); i++ )
  {
    wxArrayString ar = MyUtils::SplitString( m_points[i], "," );
    if ( ar[0] == group_name && ar[1] == point_name )
    {
      for ( int j = 0; j < 3; j++ )
        ar[j+2] = ( wxString("") << ras[j] );
      m_points[i] = MyUtils::JoinStrings( ar, "," );
      bFound = true;
      break;
    }
  }

  if ( !bFound )
  {
    wxArrayString ar;
    ar.push_back( group_name );
    ar.push_back( point_name );
    ar.push_back( ( wxString("") << ras[0] ) );
    ar.push_back( ( wxString("") << ras[1] ) );
    ar.push_back( ( wxString("") << ras[2] ) );
    m_points.push_back( MyUtils::JoinStrings( ar, "," ) );
    UpdateNames();
  }
  
  MainWindow::GetMainWindowPointer()->UpdateGotoPoints();
  
  Close();
}

void DialogSavePoint::OnDelete( wxCommandEvent& event )
{
  if ( !ValidateInput() )
    return;
  
  wxString group_name = m_comboGroupName->GetValue().Trim( true ).Trim( false );
  wxString point_name = m_comboPointName->GetValue().Trim( true ).Trim( false );
  for ( size_t i = 0; i < m_points.size(); i++ )
  {
    wxArrayString ar = MyUtils::SplitString( m_points[i], "," );
    if ( ar[0] == group_name && ar[1] == point_name )
    {
      m_points.RemoveAt( i );
      UpdateNames();
      m_comboPointName->SetValue( "" );
      break;
    }
  }
  MainWindow::GetMainWindowPointer()->UpdateGotoPoints();
}

bool DialogSavePoint::ValidateInput()
{
  wxString group_name = m_comboGroupName->GetValue().Trim( true ).Trim( false );
  wxString point_name = m_comboPointName->GetValue().Trim( true ).Trim( false );
  if ( group_name.IsEmpty() || point_name.IsEmpty() )
  {
    wxMessageDialog dlg( this, 
                         _("Group or point name can not be empty."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return false;
  }
  else
    return true;
}


void DialogSavePoint::SetGotoPoints( const wxArrayString& strgs)
{
  m_points = strgs;
  UpdateNames();
}

void DialogSavePoint::UpdateNames()
{
  wxString group_name = m_comboGroupName->GetValue().Trim( true ).Trim( false );
  m_comboGroupName->Clear();
  m_comboPointName->Clear();
  for ( size_t i = 0; i < m_points.size(); i++ )
  {
    wxArrayString ar = MyUtils::SplitString( m_points[i], "," );
    if ( m_comboGroupName->FindString( ar[0] ) == wxNOT_FOUND )
      m_comboGroupName->Append( ar[0] );
    if ( group_name == ar[0] )
      m_comboPointName->Append( ar[1] );
  }
}

void DialogSavePoint::OnComboGroupName( wxCommandEvent& event )
{
  UpdateNames();
  m_comboPointName->SetSelection( 0 );
}

void DialogSavePoint::OnShow( wxShowEvent& event )
{
//#if wxCHECK_VERSION(2,9,0)
#if wxVERSION_NUMBER > 2900
  if ( event.IsShown() )
#else
  if ( event.GetShow() )
#endif
    m_comboGroupName->SetFocus(); 
}

