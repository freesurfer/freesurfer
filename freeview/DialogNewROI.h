/**
 * @file  DialogNewROI.h
 * @brief Dialog to create new ROI.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/06/17 20:41:17 $
 *    $Revision: 1.6 $
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

