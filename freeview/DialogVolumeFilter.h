/**
 * @file  DialogVolumeFilter.h
 * @brief Dialog to create gradient volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/26 19:04:05 $
 *    $Revision: 1.1 $
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

