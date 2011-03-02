/**
 * @file  DialogVolumeFilter.h
 * @brief Dialog to create gradient volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.2 $
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
#ifndef DialogVolumeFilter_h
#define DialogVolumeFilter_h

#include <wx/wx.h>

class wxTextCtrl;
class wxSpinCtrl;
class LayerMRI;
class VolumeFilter;

class DialogVolumeFilter : public wxDialog
{
public:
  DialogVolumeFilter( wxWindow* parent );
  virtual ~DialogVolumeFilter();

  void SetFilter( VolumeFilter* filter );
  
  void OnApply( wxCommandEvent& event );
  void OnClose( wxCommandEvent& event )
  {
    Close();
  }
  
  int GetKernelSize();
  
  double GetSigma();
  
  void SetSigma( double dvalue );
  
  void ShowSigma( bool bShow );

protected:

  wxSpinCtrl*   m_spinKernelSize;
  wxButton*     m_btnApply;
  wxTextCtrl*   m_textSigma;
  
  VolumeFilter* m_filter;

  DECLARE_EVENT_TABLE()
};

#endif

