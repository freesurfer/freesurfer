/**
 * @file  DialogNewWayPoints.h
 * @brief Dialog to create new way points.
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
#ifndef DialogNewWayPoints_h
#define DialogNewWayPoints_h

#include <wx/wx.h>


class wxCheckBox;
class wxChoice;
class wxTextCtrl;
class wxRadioButton;
class LayerMRI;
class LayerCollection;

class DialogNewWayPoints : public wxDialog
{
public:
  DialogNewWayPoints( wxWindow* parent, LayerCollection* col );
  virtual ~DialogNewWayPoints();

  wxString GetWayPointsName();
  void SetWayPointsName( const wxString& name );

  LayerMRI* GetTemplate();
  
  int GetType();

  void OnOK( wxCommandEvent& event );

private:
  wxChoice*       m_choiceTemplate;
  wxTextCtrl*     m_textName;
  wxRadioButton*  m_radioWayPoints;
  wxRadioButton*  m_radioControlPoints;

  DECLARE_EVENT_TABLE()
};

#endif

