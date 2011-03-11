/**
 * @file  WorkerThread.cpp
 * @brief View for rendering.
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

#include <wx/wx.h>
#include "WorkerThread.h"
#include "LayerMRI.h"
#include "LayerROI.h"
#include "LayerDTI.h"
#include "LayerPLabel.h"
#include "LayerSurface.h"
#include "LayerWayPoints.h"
#include "LayerCollection.h"
#include "MainWindow.h"

WorkerThread::WorkerThread( wxWindow* wnd) : wxThread(),
    m_wndMain( wnd ),
    m_mri( NULL ),
    m_roi( NULL ),
    m_surface( NULL ),
    m_waypoints( NULL ),
    m_nTask( 0 ),
    m_bAllVolumes( false )
{}

bool WorkerThread::LoadVolume( LayerMRI* mri )
{
  m_mri = mri;
  m_nTask = TT_LoadVolume;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::SaveVolume( LayerMRI* mri )
{
  m_mri = mri;
  m_nTask = TT_SaveVolume;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::SaveROI( LayerROI* roi )
{
  m_roi = roi;
  m_nTask = TT_SaveROI;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::LoadSurface( LayerSurface* surface )
{
  m_surface = surface;
  m_nTask = TT_LoadSurface;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::SaveSurface( LayerSurface* surface )
{
  m_surface = surface;
  m_nTask = TT_SaveSurface;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::LoadSurfaceVector( LayerSurface* surface )
{
  m_surface = surface;
  m_nTask = TT_LoadSurfaceVector;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

bool WorkerThread::SaveWayPoints( LayerWayPoints* wp )
{
  m_waypoints = wp;
  m_nTask = TT_SaveWayPoints;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}


bool WorkerThread::RotateVolume( std::vector<RotationElement>& rotations, bool bAll )
{
  m_rotations = rotations;
  m_nTask = TT_RotateVolume;
  m_bAllVolumes = bAll;
  if ( Create() == wxTHREAD_NO_ERROR )
  {
    Run();
    return true;
  }
  else
    return false;
}

void WorkerThread::OnExit()
{}

void *WorkerThread::Entry()
{
  if ( !m_wndMain )
    return NULL;

  wxCommandEvent event( wxEVT_COMMAND_MENU_SELECTED, ID_WORKER_THREAD );
  event.SetInt( 0 );

  switch ( m_nTask )
  {
  case TT_LoadVolume:
    event.SetString( _("Load") );
    event.SetClientData( (void*)m_mri );
    break;
  case TT_SaveVolume:
    event.SetString( _("Save") );
    event.SetClientData( (void*)m_mri );
    break;
  case TT_SaveROI:
    event.SetString( _("Save") );
    event.SetClientData( (void*)m_roi );
    break;
  case TT_SaveWayPoints:
    event.SetString( _("Save") );
    event.SetClientData( (void*)m_waypoints );
    break;
  case TT_LoadSurface:
    event.SetString( _("Load") );
    event.SetClientData( (void*)m_surface );
    break;
  case TT_SaveSurface:
    event.SetString( _("Save") );
    event.SetClientData( (void*)m_surface );
    break;
  case TT_LoadSurfaceVector:
    event.SetString( _("LoadVector") );
    event.SetClientData( (void*)m_surface );
    break;
  case TT_RotateVolume:
    event.SetString( _("Rotate") );
    event.SetClientData( (void*)m_mri );
    break;
  }
  wxPostEvent( m_wndMain, event );

  bool bSuccess = false;
  try {
    switch ( m_nTask )
    {
    case TT_LoadVolume:
      if ( m_mri->IsTypeOf( "DTI" ) )
        bSuccess = ( (LayerDTI*)m_mri )->LoadDTIFromFile( m_wndMain, event );
      else if ( m_mri->IsTypeOf( "PLabel" ) )
        bSuccess = ( (LayerPLabel*)m_mri )->LoadVolumeFiles( m_wndMain, event );
      else
        bSuccess = m_mri->LoadVolumeFromFile( m_wndMain, event );
      break;
    case TT_SaveVolume:
      bSuccess = m_mri->SaveVolume( m_wndMain, event );
      break;
    case TT_LoadSurface:
      bSuccess = m_surface->LoadSurfaceFromFile( m_wndMain, event );
      break; 
    case TT_SaveSurface:
      bSuccess = m_surface->SaveSurface( m_wndMain, event );
      break;
    case TT_LoadSurfaceVector:
      bSuccess = m_surface->LoadVectorFromFile( m_wndMain, event );
      break;
    case TT_SaveROI:
      bSuccess = m_roi->SaveROI( m_wndMain, event );
      break;
    case TT_SaveWayPoints:
      bSuccess = m_waypoints->Save();
      break;
    case TT_RotateVolume:
    {
      bSuccess = true;
      std::vector<Layer*> layers = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetAllLayers();
      // first update ROI and waypoints before their reference volume is rotated
      for ( size_t i = 0; i < layers.size(); i++ )
      {
        if ( layers[i]->IsTypeOf( "ROI" ) )
        {
          ( (LayerROI*)layers[i] )->UpdateLabelData( m_wndMain, event );
        }
        else if ( layers[i]->IsTypeOf( "WayPoints" ) )
        {
          ( (LayerWayPoints*)layers[i] )->UpdateLabelData();
        }
  
      }
      
      if ( m_bAllVolumes )
      {
          // then rotate MRI volumes
        for ( size_t i = 0; i < layers.size(); i++ )
        {
          if ( layers[i]->IsTypeOf( "MRI" ) && !layers[i]->Rotate( m_rotations, m_wndMain, event ) )
          {
            bSuccess = false;
            break;
          }
        }
        // at last rotate others
        for ( size_t i = 0; i < layers.size() && bSuccess; i++ )
        {
          if ( !layers[i]->IsTypeOf( "MRI" ) && !layers[i]->Rotate( m_rotations, m_wndMain, event ) )
          {
            bSuccess = false;
            break;
          }
        }
      }
      else
      {
        LayerMRI* layer = (LayerMRI*) MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
        if ( !layer->Rotate( m_rotations, m_wndMain, event ) )
        {
          bSuccess = false;
          break;
        }
      }
      
      if ( event.GetInt() >= 100 )
        event.SetInt( 50 );
    }
    break;
    }
  }
  catch ( int nErrorCode )
  {
    wxString msg = _("Failed");
    if ( event.GetString() == _("Load") )
    {
      Layer* layer = (Layer* )(void* )event.GetClientData();
      delete layer;
    }
    event.SetString( msg );
    wxPostEvent( m_wndMain, event );
    return NULL;
  }

  if ( bSuccess )
  {
    event.SetInt( -1 );
  }
  else
  {
    wxString msg = _("Failed");
    if ( m_mri )
      msg += wxString::FromAscii( m_mri->GetErrorString().c_str() );

    if ( event.GetString() == _("Load") )
    {
      Layer* layer = (Layer* )(void* )event.GetClientData();
      delete layer;
    }

    event.SetString( msg );
  }
  
  wxPostEvent( m_wndMain, event );

  return NULL;
}

