/**
 * @file  DialogNewROI.h
 * @brief Dialog to create new ROI.
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
#ifndef DialogNewROI_h
#define DialogNewROI_h

#include <wx/wx.h>


class wxCheckBox;
class wxChoice;
class wxTextCtrl;
class LayerMRI;
class LayerCollection;

class DialogNewROI : public wxDialog
{
public:
  DialogNewROI( wxWindow* parent, LayerCollection* col );
  virtual ~DialogNewROI();

  wxString GetROIName();
  void SetROIName( const wxString& name );

  LayerMRI* GetTemplate();

  void OnOK( wxCommandEvent& event );

  void OnTextEnter( wxCommandEvent& event );

private:
  wxChoice*  m_choiceTemplate;
  wxTextCtrl*  m_textName;

  DECLARE_EVENT_TABLE()
};

#endif

