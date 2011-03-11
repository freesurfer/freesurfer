/**
 * @file  BuildContourThread.h
 * @brief Worker thread class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:35 $
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

#ifndef BuildContourThread_h
#define BuildContourThread_h


#include <wx/thread.h>

class wxWindow;
class LayerMRI;

class BuildContourThread : public wxThread
{
public:
  BuildContourThread( wxWindow* wnd );

  bool BuildContour( LayerMRI* mri, int nSegValue, int nThreadID );
  virtual void* Entry();

  virtual void OnExit();
  
  void Abort();
  
  void SetSmoothFactor( int nFactor )
  {
    m_nSmoothFactor = nFactor;
  }

protected:
  wxWindow*   m_wndMain;
  LayerMRI*   m_mri;
  bool        m_bAbort;
  int         m_nSegValue;
  int         m_nThreadID;
  int         m_nSmoothFactor;
};

#endif


