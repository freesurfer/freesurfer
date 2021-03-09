/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "ThreadIOWorker.h"
#include "LayerMRI.h"
#include "LayerSurface.h"
#include "LayerDTI.h"
#include "LayerPLabel.h"
#include "LayerTrack.h"
#include "LayerVolumeTrack.h"
#include "LayerConnectomeMatrix.h"
#include "LayerFCD.h"
#include "LayerODF.h"
#include <QApplication>
#include <QDebug>

ThreadIOWorker::ThreadIOWorker(QObject *parent) :
  QThread(parent),
  m_nJobType( JT_LoadVolume ),
  m_layer( NULL )
{
}

void ThreadIOWorker::LoadVolume( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_LoadVolume;
  m_args = args;
  start();
}

void ThreadIOWorker::SaveVolume( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_SaveVolume;
  m_args = args;
  start();
}

void ThreadIOWorker::TransformVolume( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_TransformVolume;
  m_args = args;
  start();
}

void ThreadIOWorker::LoadSurface( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_LoadSurface;
  m_args = args;
  start();
}

void ThreadIOWorker::SaveSurface( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_SaveSurface;
  m_args = args;
  start();
}


void ThreadIOWorker::LoadSurfaceOverlay( Layer* layer, const QVariantMap& args )
{
  m_layer = layer;
  m_nJobType = JT_LoadSurfaceOverlay;
  m_args = args;
  start();
}

void ThreadIOWorker::LoadTrack( Layer* layer, const QVariantMap& args)
{
  m_layer = layer;
  m_nJobType = JT_LoadTrack;
  m_args = args;
  start();
}

void ThreadIOWorker::LoadConnectomeMatrix(Layer* layer, const QVariantMap& args)
{
  m_layer = layer;
  m_nJobType = JT_LoadConnectome;
  m_args = args;
  start();
}

void ThreadIOWorker::LoadFCD(Layer* layer, const QVariantMap& args)
{
  m_layer = layer;
  m_nJobType = JT_LoadFCD;
  m_args = args;
  start();
}

void ThreadIOWorker::LoadODF(Layer *layer, const QVariantMap &args)
{
  m_layer = layer;
  m_nJobType = JT_LoadODF;
  m_args = args;
  start();
}

void ThreadIOWorker::run()
{
  connect(qApp, SIGNAL(GlobalProgress(int)), this, SIGNAL(Progress(int)), Qt::UniqueConnection);
  if ( m_nJobType == JT_LoadVolume )
  {
    if (m_layer->IsTypeOf("DTI"))
    {
      LayerDTI* mri = qobject_cast<LayerDTI*>( m_layer );
      if ( !mri )
      {
        return;
      }
      if ( !mri->LoadDTIFromFile() )
      {
        emit Error( m_layer, m_nJobType );
      }
      else
      {
        emit Finished( m_layer, m_nJobType );
      }
    }
    else if (m_layer->IsTypeOf("PLabel"))
    {
      LayerPLabel* mri = qobject_cast<LayerPLabel*>( m_layer );
      if ( !mri )
      {
        return;
      }
      if ( !mri->LoadVolumeFiles() )
      {
        emit Error( m_layer, m_nJobType );
      }
      else
      {
        emit Finished( m_layer, m_nJobType );
      }
    }
    else if (m_layer->IsTypeOf("VolumeTrack"))
    {
      LayerVolumeTrack* mri = qobject_cast<LayerVolumeTrack*>( m_layer );
      if ( !mri )
      {
        return;
      }
      if ( !mri->LoadFromFile() )
      {
        emit Error( m_layer, m_nJobType );
      }
      else
      {
        emit Finished( m_layer, m_nJobType );
      }
    }
    else
    {
      LayerMRI* mri = qobject_cast<LayerMRI*>( m_layer );
      if ( !mri )
      {
        return;
      }
      if ( !mri->LoadVolumeFromFile() )
      {
        emit Error( m_layer, m_nJobType );
      }
      else
      {
        emit Finished( m_layer, m_nJobType );
      }
    }
  }
  else if (m_nJobType == JT_SaveVolume)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>( m_layer );
    if ( !mri )
    {
      return;
    }
    if ( !mri->SaveVolume() )
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      emit Finished( m_layer, m_nJobType );
    }
  }
  else if (m_nJobType == JT_TransformVolume)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>( m_layer );
    if ( !mri )
    {
      return;
    }
    if (m_args.value("unload").toBool())
    {
      mri->UnloadVolumeTransform();
      emit Finished( m_layer, m_nJobType );
    }
    else if ( !mri->LoadVolumeTransform() )
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      emit Finished( m_layer, m_nJobType );
    }
  }
  else if ( m_nJobType == JT_LoadSurface )
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>( m_layer );
    if ( !surf )
    {
      return;
    }
    if ( !surf->LoadSurfaceFromFile(m_args["ignore_vg"].toBool()) )
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      if (m_args.value("hidden").toBool())
        m_layer->setProperty("hidden", true);

      emit Finished( m_layer, m_nJobType );
    }
  }
  else if ( m_nJobType == JT_SaveSurface )
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>( m_layer );
    if ( !surf )
    {
      return;
    }
    if ( !surf->SaveSurface() )
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      emit Finished( m_layer, m_nJobType );
    }
  }
  else if ( m_nJobType == JT_LoadSurfaceOverlay )
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>( m_layer );
    if ( !surf )
    {
      return;
    }
    QString fn = m_args["FileName"].toString();
    QString fn_reg = m_args["Registration"].toString();
    bool bCorrelation = m_args["Correlation"].toBool();
    bool bSecondHalf = m_args["SecondHalfData"].toBool();
    if ( !surf->LoadOverlayFromFile(fn, fn_reg, bCorrelation, bSecondHalf))
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      emit Finished( m_layer, m_nJobType );
    }
  }
  else if ( m_nJobType == JT_LoadTrack )
  {
    LayerTrack* layer = qobject_cast<LayerTrack*>( m_layer );
    if ( !layer )
    {
      return;
    }
    connect(layer, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)), Qt::UniqueConnection);
    if ( !layer->LoadTrackFromFiles() )
    {
      emit Error( m_layer, m_nJobType );
    }
    else
    {
      emit Finished( m_layer, m_nJobType );
    }
  }
  else if (m_nJobType == JT_LoadConnectome)
  {
    LayerConnectomeMatrix* layer = qobject_cast<LayerConnectomeMatrix*>(m_layer);
    if (!layer)
      return;
    if (!layer->LoadFromFile())
    {
      emit Error(m_layer, m_nJobType);
    }
    else
    {
      emit Finished(m_layer, m_nJobType);
    }
  }
  else if (m_nJobType == JT_LoadFCD)
  {
    LayerFCD* layer = qobject_cast<LayerFCD*>(m_layer);
    if (!layer)
      return;
    if (!layer->Load(m_args["SubjectDir"].toString(), m_args["Subject"].toString(), m_args["Suffix"].toString()))
    {
      emit Error(m_layer, m_nJobType);
    }
    else
    {
      emit FCDLoadFinished(layer);
    }
  }
  else if (m_nJobType == JT_LoadODF)
  {
    LayerODF* layer = qobject_cast<LayerODF*>(m_layer);
    if (!layer)
      return;
    if (!layer->Load(m_args["Filename"].toString()))
    {
      emit Error(m_layer, m_nJobType);
    }
    else
    {
      emit Finished(m_layer, m_nJobType);
    }
  }
  disconnect(qApp, SIGNAL(GlobalProgress(int)), this, SIGNAL(Progress(int)));
}
