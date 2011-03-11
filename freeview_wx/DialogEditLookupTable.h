/**
 * @file  DialogEditLookupTable.h
 * @brief Dialog window for editing lookup table.
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
#ifndef DialogEditLookupTable_h
#define DialogEditLookupTable_h

#include <wx/wx.h>
#include "LUTDataHolder.h"

class wxTextCtrl;
class wxChoice;
class wxListBox;
class wxColourPickerCtrl;
class wxComboBox;
class wxColorIndicator;
class LUTDataHolder;

class DialogEditLookupTable : public wxDialog
{
public:
  DialogEditLookupTable( wxWindow* parent );
  virtual ~DialogEditLookupTable();

protected:
  void OnButtonClose  ( wxCommandEvent& event );
  void OnButtonAddItem ( wxCommandEvent& event );
  void OnButtonUpdateItem ( wxCommandEvent& event );
  void OnButtonRemoveItem ( wxCommandEvent& event );
  void OnButtonImportItem ( wxCommandEvent& event );
  void OnButtonNewTable ( wxCommandEvent& event );
  void OnButtonSaveTable ( wxCommandEvent& event );
  void OnButtonDeleteTable( wxCommandEvent& event );

  void OnChoiceExistingTables ( wxCommandEvent& event );
  void OnComboEditableTables ( wxCommandEvent& event );

  void OnListExistingTable ( wxCommandEvent& event );
  void OnListNewTable   ( wxCommandEvent& event );

private:
  void Initialize();
  void PopulateColorTable( COLOR_TABLE* ct, wxListBox* listbox );

  wxButton*  m_btnAddItem;
  wxButton*  m_btnUpdateItem;
  wxButton*  m_btnRemoveItem;
  wxButton*  m_btnImportItem;
  wxButton*  m_btnNewTable;
  wxButton*  m_btnSaveTable;
  wxButton*  m_btnDeleteTable;
  wxColorIndicator* m_colorIndicator;
  wxColourPickerCtrl* m_colorPicker;
  wxTextCtrl*  m_textLabel;
  wxComboBox*  m_comboEditableTables;
  wxChoice*  m_choiceExistingTables;
// wxTextCtrl*  m_textFA;
  wxListBox*  m_listExistingTable;
  wxListBox*  m_listNewTable;

  wxString  m_strLastDir;

  LUTDataHolder* m_luts;
  COLOR_TABLE* m_ctTemplate;
  COLOR_TABLE* m_ctEdit;

  DECLARE_EVENT_TABLE()
};

#endif

