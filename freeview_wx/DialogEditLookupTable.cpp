/**
 * @file  DialogEditLookupTable.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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



#include "DialogEditLookupTable.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/clrpicker.h>
#include "MainWindow.h"
#include "wxColorIndicator.h"

BEGIN_EVENT_TABLE( DialogEditLookupTable, wxDialog )
EVT_BUTTON ( XRCID( "ID_BUTTON_CLOSE" ),  DialogEditLookupTable::OnButtonClose )
EVT_BUTTON ( XRCID( "ID_BUTTON_ADD_ITEM" ),  DialogEditLookupTable::OnButtonAddItem )
EVT_BUTTON ( XRCID( "ID_BUTTON_UPDATE_ITEM" ), DialogEditLookupTable::OnButtonUpdateItem )
EVT_BUTTON ( XRCID( "ID_BUTTON_REMOVE_ITEM" ), DialogEditLookupTable::OnButtonRemoveItem )
EVT_BUTTON ( XRCID( "ID_BUTTON_IMPORT_ITEM" ), DialogEditLookupTable::OnButtonImportItem )
EVT_BUTTON ( XRCID( "ID_BUTTON_NEW_TABLE" ), DialogEditLookupTable::OnButtonNewTable )
EVT_BUTTON ( XRCID( "ID_BUTTON_SAVE_TABLE" ), DialogEditLookupTable::OnButtonSaveTable )
EVT_BUTTON ( XRCID( "ID_BUTTON_DELETE_TABLE" ), DialogEditLookupTable::OnButtonDeleteTable )
EVT_CHOICE ( XRCID( "ID_CHOICE_EXISTING_TABLES" ), DialogEditLookupTable::OnChoiceExistingTables )
EVT_COMBOBOX ( XRCID( "ID_COMBO_EDITABLE_TABLES" ),  DialogEditLookupTable::OnComboEditableTables )
EVT_LISTBOX ( XRCID( "ID_LISTBOX_EXISTING_TABLE" ),  DialogEditLookupTable::OnListExistingTable )
EVT_LISTBOX ( XRCID( "ID_LISTBOX_NEW_TABLE" ),  DialogEditLookupTable::OnListNewTable )
END_EVENT_TABLE()


DialogEditLookupTable::DialogEditLookupTable( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_EDIT_LUT") );
  m_choiceExistingTables = 
    XRCCTRL( *this, "ID_CHOICE_EXISTING_TABLES", wxChoice );
  m_listExistingTable = 
    XRCCTRL( *this, "ID_LISTBOX_EXISTING_TABLE", wxListBox );
  m_listNewTable = XRCCTRL( *this, "ID_LISTBOX_NEW_TABLE", wxListBox );
  m_colorIndicator = XRCCTRL( *this, "ID_COLORINDICATOR", wxColorIndicator );
  m_colorPicker = XRCCTRL( *this, "ID_COLORPICKER", wxColourPickerCtrl );
  m_textLabel = XRCCTRL( *this, "ID_TEXT_LABEL", wxTextCtrl );

  Initialize();
}

DialogEditLookupTable::~DialogEditLookupTable()
{}

void DialogEditLookupTable::Initialize()
{
  m_luts = MainWindow::GetMainWindowPointer()->GetLUTData();
  for ( int i = 0; i < m_luts->GetCount(); i++ )
  {
    m_choiceExistingTables->Append
      ( wxString::FromAscii(m_luts->GetName( i )), 
	(void*)m_luts->GetColorTable( i ) );
    if ( i == 0 )
    {
      m_ctTemplate = m_luts->GetColorTable( 0 );
      PopulateColorTable( m_ctTemplate, m_listExistingTable );
    }
  }

  m_choiceExistingTables->SetSelection( 0 );
}

void DialogEditLookupTable::OnButtonClose( wxCommandEvent& event )
{
  /* if ( GetVectorFileName().IsEmpty() )
   {
    wxMessageDialog dlg
    ( this, 
    _("Vector file name can not be empty."), 
    _("Error"), 
    wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
   }
   else if ( GetFAFileName().IsEmpty() )
   {
    wxMessageDialog dlg
    ( this, 
    _("FA file name can not be empty."), 
    _("Error"), 
    wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }*/

  Close();
}

void DialogEditLookupTable::PopulateColorTable( COLOR_TABLE* ct, 
						wxListBox* listbox )
{
  listbox->Clear();
  listbox->FitInside();
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( ct, &nTotalCount );
  int nValid;
  char name[1000];
  int nValidCount = 0;
  for ( int i = 0; i < nTotalCount; i++ )
  {
    CTABisEntryValid( ct, i, &nValid );
    if ( nValid )
    {
      CTABcopyName( ct, i, name, 1000 );
      listbox->Append( wxString::Format( _("%d: %s"), i, name ) );

      nValidCount++;
    }
  }
}

void DialogEditLookupTable::OnListExistingTable( wxCommandEvent& event )
{
  if ( m_ctTemplate )
  {
    wxString strg = m_listExistingTable->GetStringSelection();
    strg = strg.Left( strg.Find( ':' ) );
    long nIndex;
    int nr, ng, nb;
    if ( strg.ToLong( &nIndex ) && 
	 CTABrgbAtIndexi( m_ctTemplate, nIndex, &nr, &ng, &nb ) == 0 )
    {
      m_colorIndicator->SetColor( wxColour( nr, ng, nb ) );
    }
  }
}

void DialogEditLookupTable::OnListNewTable( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonImportItem( wxCommandEvent& event )
{
  wxArrayInt sels;
  m_listExistingTable->GetSelections( sels );
  for ( int i = 0; i < (int)sels.GetCount(); i++ )
  {
    // wxString strg = m_listExistingTable->GetSelectionString( sels[i] );

  }
}

void DialogEditLookupTable::OnButtonAddItem( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonUpdateItem( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonRemoveItem( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonNewTable( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonSaveTable( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnButtonDeleteTable( wxCommandEvent& event )
{}

void DialogEditLookupTable::OnChoiceExistingTables( wxCommandEvent& event )
{
  COLOR_TABLE* ct = ( COLOR_TABLE* )( void* )m_choiceExistingTables->
    GetClientData( event.GetSelection() );
  if ( ct )
  {
    m_ctTemplate = ct;
    PopulateColorTable( ct, m_listExistingTable );
  }
}

void DialogEditLookupTable::OnComboEditableTables( wxCommandEvent& event )
{}

