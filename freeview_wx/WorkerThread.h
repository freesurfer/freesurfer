/**
 * @file  WorkerThread.h
 * @brief Worker thread class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
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

#ifndef WorkerThread_h
#define WorkerThread_h


#include <wx/thread.h>
#include "CommonDataStruct.h"

class wxWindow;
class LayerMRI;
class LayerROI;
class LayerDTI;
class LayerSurface;
class LayerWayPoints;

class WorkerThread : public wxThread
{
public:
  WorkerThread( wxWindow* wnd );

  bool LoadVolume( LayerMRI* mri );
  bool SaveVolume( LayerMRI* mri );
  bool SaveROI( LayerROI* roi );
  bool LoadSurface( LayerSurface* surface );
  bool LoadSurfaceVector( LayerSurface* surface );
  bool SaveSurface( LayerSurface* surface );
  bool SaveWayPoints( LayerWayPoints* wp );
  bool RotateVolume( std::vector<RotationElement>& rotations, bool bAllVolumes );

  virtual void* Entry();

  virtual void OnExit();

protected:
  enum TaskType { TT_LoadVolume = 0,
                  TT_SaveVolume,
                  TT_LoadSurface,
                  TT_SaveSurface,
                  TT_LoadSurfaceVector,
                  TT_SaveROI,
                  TT_SaveWayPoints,
                  TT_RotateVolume };

  wxWindow* m_wndMain;
  LayerMRI* m_mri;
  LayerROI* m_roi;
  LayerSurface* m_surface;
  LayerWayPoints* m_waypoints;

  std::vector<RotationElement> m_rotations;

  int   m_nTask;
  bool  m_bAllVolumes;
};

#endif


