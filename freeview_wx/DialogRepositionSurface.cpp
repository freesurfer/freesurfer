/**
 * @file  DialogRepositionSurface.h
 * @brief Dialog to create gradient volume.
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



#include "DialogRepositionSurface.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "LayerMRI.h"
#include "LayerSurface.h"
#include "MainWindow.h"
#include "MyUtils.h"

BEGIN_EVENT_TABLE( DialogRepositionSurface, wxDialog )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_APPLY" ),     DialogRepositionSurface::OnApply )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CLOSE" ),     DialogRepositionSurface::OnClose )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_UNDO" ),      DialogRepositionSurface::OnUndo )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE" ),      DialogRepositionSurface::OnSave )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_SAVE_AS" ),   DialogRepositionSurface::OnSaveAs )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_TARGET" ),    DialogRepositionSurface::OnChoiceTarget )
END_EVENT_TABLE()


DialogRepositionSurface::DialogRepositionSurface( wxWindow* parent ) : Listener( "DialogRepositionSurface" )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_REPOSITION_SURFACE") );
  m_textVertex    =   XRCCTRL( *this, "ID_TEXT_VERTEX",     wxTextCtrl );
  m_textSize      =   XRCCTRL( *this, "ID_TEXT_SIZE",       wxTextCtrl );
  m_textTarget    =   XRCCTRL( *this, "ID_TEXT_TARGET",     wxTextCtrl );
  m_textSigma     =   XRCCTRL( *this, "ID_TEXT_SIGMA",      wxTextCtrl );
  m_choiceTarget  =   XRCCTRL( *this, "ID_CHOICE_TARGET",   wxChoice );
  m_btnSave       =   XRCCTRL( *this, "ID_BUTTON_SAVE",     wxButton );
  m_btnSaveAs     =   XRCCTRL( *this, "ID_BUTTON_SAVE_AS",  wxButton );
  m_btnUndo       =   XRCCTRL( *this, "ID_BUTTON_UNDO",     wxButton );
}

DialogRepositionSurface::~DialogRepositionSurface()
{
}

void DialogRepositionSurface::OnApply( wxCommandEvent& event )
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "Surface" );
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
  wxString msg;
  if ( !surf )
    msg = _("No active surface found.");
  else if ( !mri )
    msg = _("No active volume found." );
  
  if ( !msg.IsEmpty() )
  {
    wxMessageDialog dlg( this, msg, _("Error"), wxOK );
    dlg.ShowModal();
    return;
  }
  if ( ValidateAll() )
  {
    if ( m_choiceTarget->GetCurrentSelection() == 0 )
    {
      surf->RepositionSurface( mri, GetVertex(), 
                           GetIntensity(),
                           GetNeighborSize(),
                           GetSigma() );
    }
    else
    {
      double pos[3];
      GetCoordinate( pos );
      surf->RepositionSurface( mri, GetVertex(), 
                               pos,
                               GetNeighborSize(),
                               GetSigma() );
    }
    UpdateUI();
  }
}

void DialogRepositionSurface::OnChoiceTarget( wxCommandEvent& event )
{
  m_textTarget->SetSize( (event.GetSelection() == 0 ? 80 : 140), m_textTarget->GetSize().GetHeight() );
}

void DialogRepositionSurface::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerModified" )
  {
    UpdateUI();
  }
  else if ( iMsg == "SurfaceVertexClicked" )
  {   
    LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "Surface" );
    if ( surf )
    {
      int nVertex = surf->GetVertexIndexAtTarget( surf->GetSlicePosition(), NULL );
      if ( nVertex >= 0 )
        m_textVertex->SetValue( ( wxString() << nVertex ) );
    }
  }
}

void DialogRepositionSurface::UpdateUI()
{
  MainWindow* mainwnd = MainWindow::GetMainWindowPointer();
  LayerSurface* surf = (LayerSurface*)mainwnd->GetActiveLayer( "Surface" );
  m_btnSave->Enable( surf && surf->IsModified() && !mainwnd->IsProcessing() );
  m_btnSaveAs->Enable( surf && !mainwnd->IsProcessing() );
  m_btnUndo->Enable( surf && surf->HasUndo() );
}

void DialogRepositionSurface::OnSave( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveSurface();
  UpdateUI();
}

void DialogRepositionSurface::OnSaveAs( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SaveSurfaceAs();
  UpdateUI();
}

int DialogRepositionSurface::GetVertex()
{
  long value;
  m_textVertex->GetValue().ToLong( &value );
  return value;
}

int DialogRepositionSurface::GetNeighborSize()
{
  long value;
  m_textSize->GetValue().ToLong( &value );
  return value;
}

double DialogRepositionSurface::GetIntensity()
{
  double value;
  m_textTarget->GetValue().ToDouble( &value );
  return value;
}

void DialogRepositionSurface::GetCoordinate( double* pos )
{
  wxArrayString sa = MyUtils::SplitString( m_textTarget->GetValue(), "," );
  if ( sa.size() < 3 )
    sa = MyUtils::SplitString( m_textTarget->GetValue(), " " );
  for ( int i = 0; i < 3; i++ )
    sa[i].ToDouble( pos+i );
}

double DialogRepositionSurface::GetSigma()
{
  double value;
  m_textSigma->GetValue().ToDouble( &value );
  return value;
}

void DialogRepositionSurface::OnUndo( wxCommandEvent& event )
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "Surface" );
  if ( surf )
  {
    surf->Undo();
    UpdateUI();
  }
}

bool DialogRepositionSurface::ValidateAll()
{
  wxString name;
  long nval;
  double dval;
  if ( !m_textVertex->GetValue().ToLong( &nval ) )
    name = _("Vertex");
  else if ( !m_textSize->GetValue().ToLong( &nval ) )
    name = _("Size");
  else if ( !m_textSigma->GetValue().ToDouble( &dval ) )
    name = _("Sigma");
  else if ( m_choiceTarget->GetCurrentSelection() == 0 && !m_textTarget->GetValue().ToDouble( &dval ) )
    name = _("Intensity");
  else if ( m_choiceTarget->GetCurrentSelection() == 1 )
  {
    wxArrayString sa = MyUtils::SplitString( m_textTarget->GetValue(), "," );
    if ( sa.size() < 3 )
      sa = MyUtils::SplitString( m_textTarget->GetValue(), " " );
    if ( sa.size() < 3 || !sa[0].ToDouble( &dval ) || !sa[1].ToDouble( &dval ) || !sa[2].ToDouble( &dval ) )
      name = _("Coordinate");
  }
  
  if ( !name.IsEmpty() )
  {
    wxMessageDialog dlg( this, wxString("Invalid input for ") + name, 
                        _("Error"), wxOK );
    dlg.ShowModal();
    return false;
  }
  else
    return true;
}
