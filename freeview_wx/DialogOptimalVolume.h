/**
 * @file  DialogOptimalVolume.h
 * @brief Dialog to compute optimal volume.
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
#ifndef DialogOptimalVolume_h
#define DialogOptimalVolume_h

#include <wx/wx.h>
#include <vector>

class wxCheckBox;
class wxChoice;
class wxTextCtrl;
class wxCheckListBox;
class LayerMRI;
class LayerCollection;

class DialogOptimalVolume : public wxDialog
{
public:
  DialogOptimalVolume( wxWindow* parent, LayerCollection* col );
  virtual ~DialogOptimalVolume();

  wxString GetVolumeName();
  void SetVolumeName( const wxString& name );

  LayerMRI* GetLabelVolume();

  std::vector<LayerMRI*> GetSelectedLayers();

  void OnOK( wxCommandEvent& event );

  void OnTextEnter( wxCommandEvent& event );

private:
  wxChoice*  m_choiceLabelVolume;
  wxTextCtrl*  m_textName;
  wxCheckListBox* m_listBoxLayers;

  DECLARE_EVENT_TABLE()
};

#endif

