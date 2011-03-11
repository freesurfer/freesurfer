/**
 * @file  DialogGradientVolume.h
 * @brief Dialog to create gradient volume.
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
#ifndef DialogGradientVolume_h
#define DialogGradientVolume_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;
class LayerMRI;
class VolumeFilterGradient;

class DialogGradientVolume : public wxDialog
{
public:
  DialogGradientVolume( LayerMRI* input, LayerMRI* output, wxWindow* parent );
  virtual ~DialogGradientVolume();

  void OnApply( wxCommandEvent& event );
  void OnClose( wxCommandEvent& event )
  {
    Close();
  }
  
  void SetVolumes( LayerMRI* input, LayerMRI* output );

protected:
  void OnCheckSmooth( wxCommandEvent& event );

  wxTextCtrl*   m_textSD;
  wxCheckBox*   m_checkSmooth;
  wxButton*     m_btnApply;
  
  VolumeFilterGradient* m_filter;

  DECLARE_EVENT_TABLE()
};

#endif

