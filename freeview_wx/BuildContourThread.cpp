/**
 * @file  BuildContourThread.cpp
 * @brief View for rendering.
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

#include <wx/wx.h>
#include "BuildContourThread.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "MyUtils.h"
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"

BuildContourThread::BuildContourThread( wxWindow* wnd) : wxThread(),
    m_wndMain( wnd ),
    m_mri( NULL ),
    m_bAbort( false ),
    m_nSmoothFactor( 0 )
{}

bool BuildContourThread::BuildContour( LayerMRI* mri, int nSegValue, int nThreadID )
{
  m_mri = mri;
  m_nSegValue = nSegValue;
  m_nThreadID = nThreadID;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run(); 
    return true;
  }
  else
  {
    return false;
  }
}

void BuildContourThread::OnExit()
{}

void BuildContourThread::Abort()
{
  m_bAbort = true;
}

void *BuildContourThread::Entry()
{
  if ( !m_wndMain )
    return NULL;

  wxCommandEvent event( wxEVT_COMMAND_MENU_SELECTED, ID_THREAD_BUILD_CONTOUR );

  double dTh1 = m_mri->GetProperties()->GetContourMinThreshold();
  double dTh2 = m_mri->GetProperties()->GetContourMaxThreshold();
  bool bExtractAllRegions = m_mri->GetProperties()->GetContourExtractAllRegions();
  if ( m_nSegValue >= 0 )
  {
    dTh1 = m_nSegValue - 0.5;
    dTh2 = m_nSegValue + 0.5;
  }

  vtkActor* actor = vtkActor::New();
  actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
//  int ext[6] = { 370, 520, 0, 140, 0, 280 };
  MyUtils::BuildContourActor( m_mri->GetImageData(), dTh1, dTh2, actor, m_nSmoothFactor, NULL, bExtractAllRegions );
  m_mri->m_actorContourTemp = actor;
  actor->Delete();
  
  if ( !m_bAbort )
  {
    event.SetInt( m_nThreadID );
    event.SetClientData( (void*)m_mri );
    wxPostEvent( m_wndMain, event );
  }

  return NULL;
}

