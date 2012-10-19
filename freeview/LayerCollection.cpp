/**
 * @file  LayerCollection.cpp
 * @brief Collection of layers of the same type.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/19 15:52:08 $
 *    $Revision: 1.33 $
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
 *
 */

#include "LayerCollection.h"
#include "Layer.h"
#include "LayerProperty.h"
#include <math.h>
#include <QDebug>
#include <iostream>

LayerCollection::LayerCollection( const QString& strType, QObject* parent ) :
  QObject( parent ),
  m_layerActive( NULL ),
  m_strType( strType )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = 0;
    m_dWorldOrigin[i] = 0;
    m_dWorldVoxelSize[i] = 1;
  }
}

LayerCollection::~LayerCollection()
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    delete m_layers[i];
  }
}

bool LayerCollection::IsEmpty()
{
  return m_layers.size() == 0;
}

int LayerCollection::GetLayerIndex( Layer* layer )
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i] == layer )
    {
      return i;
    }
  }
  return -1;
}

bool LayerCollection::AddLayer( Layer* layer, bool initializeCoordinate )
{
  if ( !layer->IsTypeOf( m_strType ) )
  {
    std::cerr << "Can not add layer type of " << qPrintable(layer->GetEndType())
              << " to layer collection type of " <<  qPrintable(m_strType) << "\n";
    return false;
  }

  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i] == layer )
    {
      return false;
    }
  }

  if ( initializeCoordinate)
  {
    layer->GetSlicePosition( m_dSlicePosition );
    layer->GetWorldOrigin( m_dWorldOrigin );
    layer->GetWorldSize( m_dWorldSize );
    layer->GetWorldVoxelSize( m_dWorldVoxelSize );
  }
  else
  {
    layer->SetSlicePosition( m_dSlicePosition );
  }

  m_layers.insert( m_layers.begin(), layer );
  connect( layer, SIGNAL(ActorUpdated()), this, SIGNAL(LayerActorUpdated()) );
  connect( layer, SIGNAL(Transformed()), this, SIGNAL(LayerActorUpdated()) );
  connect( layer, SIGNAL(ActorChanged()), this, SIGNAL(LayerActorChanged()) );
  connect( layer, SIGNAL(NameChanged(QString)), this, SIGNAL(LayerNameChanged()));
  if (layer->IsTypeOf("Editable"))
    connect( layer, SIGNAL(Modified()), this, SIGNAL(LayerModified()));
  if (layer->GetProperty())
  {
    connect( layer->GetProperty(), SIGNAL(PropertyChanged()), this, SIGNAL(LayerPropertyChanged()));
    connect( layer->GetProperty(), SIGNAL(ShowInfoChanged(bool)), this, SIGNAL(LayerShowInfoChanged()));
  }
  connect( layer, SIGNAL(VisibilityChanged(bool)), this, SIGNAL(LayerVisibilityChanged()));

  this->SetActiveLayer( layer );
  emit LayerAdded( layer );

  return true;
}

bool LayerCollection::RemoveLayer( Layer* layer, bool deleteObject )
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i] == layer )
    {
      m_layers.erase( m_layers.begin() + i );

      if ( m_layers.size() == 0 )
      {
        SetActiveLayer( NULL );
      }
      else
      {
        if ( i == m_layers.size() )
        {
          SetActiveLayer( m_layers[m_layers.size()-1] );
        }
        else
        {
          SetActiveLayer( m_layers[i] );
        }
      }

      emit LayerRemoved( layer );

      if (deleteObject)
      {
        delete layer;
      }

      return true;
    }
  }

  return false;
}

void LayerCollection::MoveLayerUp()
{
  if ( this->m_layerActive )
  {
    MoveLayerUp( this->m_layerActive );
  }
}

void LayerCollection::MoveLayerDown()
{
  if ( this->m_layerActive )
  {
    MoveLayerDown( this->m_layerActive );
  }
}

bool LayerCollection::MoveLayerUp( Layer* layer )
{
  QList<Layer*> unlocked_layers;
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( !m_layers[i]->IsLocked() )
    {
      unlocked_layers << m_layers[i];
    }
  }

  for ( int i = 1; i < unlocked_layers.size(); i++)
  {
    if ( unlocked_layers[i] == layer )
    {
      Layer* temp = unlocked_layers[i-1];
      unlocked_layers[i-1] = layer;
      unlocked_layers[i] = temp;

      // restore locked layers
      for ( int j = 0; j < m_layers.size(); j++ )
      {
        if ( m_layers[j]->IsLocked() )
        {
          if ( j < unlocked_layers.size() )
          {
            unlocked_layers.insert( unlocked_layers.begin() + j, m_layers[j] );
          }
          else
          {
            unlocked_layers.push_back( m_layers[j] );
          }
        }
      }
      m_layers = unlocked_layers;

      emit LayerMoved( layer );

      return true;
    }
  }
  return false;
}

bool LayerCollection::MoveLayerDown( Layer* layer )
{
  QList<Layer*> unlocked_layers;
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( !m_layers[i]->IsLocked() )
    {
      unlocked_layers.push_back( m_layers[i] );
    }
  }

  for ( int i = 0; i < unlocked_layers.size()-1; i++)
  {
    if ( unlocked_layers[i] == layer )
    {
      Layer* temp = unlocked_layers[i+1];
      unlocked_layers[i+1] = layer;
      unlocked_layers[i] = temp;

      // restore locked layers
      for ( int j = 0; j < m_layers.size(); j++ )
      {
        if ( m_layers[j]->IsLocked() )
        {
          if ( j < unlocked_layers.size() )
          {
            unlocked_layers.insert( unlocked_layers.begin() + j, m_layers[j] );
          }
          else
          {
            unlocked_layers.push_back( m_layers[j] );
          }
        }
      }
      m_layers = unlocked_layers;

      emit LayerMoved( layer );

      return true;
    }
  }
  return false;
}

bool LayerCollection::MoveToTop( Layer* layer )
{
  QList<Layer*> unlocked_layers;
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( !m_layers[i]->IsLocked() )
    {
      unlocked_layers.push_back( m_layers[i] );
    }
  }

  for ( int i = 0; i < unlocked_layers.size(); i++)
  {
    if ( unlocked_layers[i] == layer )
    {
      for ( int j = i; j > 0 ; j-- )
      {
        unlocked_layers[j] = unlocked_layers[j-1];
      }
      unlocked_layers[0] = layer;

      // restore locked layers
      for ( int j = 0; j < m_layers.size(); j++ )
      {
        if ( m_layers[j]->IsLocked() )
        {
          if ( j < unlocked_layers.size() )
          {
            unlocked_layers.insert( unlocked_layers.begin() + j, m_layers[j] );
          }
          else
          {
            unlocked_layers.push_back( m_layers[j] );
          }
        }
      }
      m_layers = unlocked_layers;

      emit LayerMoved( layer );

      return true;
    }
  }
  return false;
}

bool LayerCollection::CycleLayer( bool bMoveUp )
{
  if ( (int)m_layers.size() > 1 )
  {
    int nActive = GetLayerIndex( m_layerActive );

    // first get unlocked layers only
    QList<Layer*> unlocked_layers;
    for ( int i = 0; i < m_layers.size(); i++ )
    {
      if ( !m_layers[i]->IsLocked() )
      {
        unlocked_layers.push_back( m_layers[i] );
      }
    }

    // record the visibilities of each layer before cycling
    bool* bVisibility = new bool[m_layers.size()];
    for ( int i = 0; i < m_layers.size(); i++ )
    {
      bVisibility[i] = m_layers[i]->IsVisible();
    }

    if ( unlocked_layers.size() == 0 )
    {
      delete[] bVisibility;
      return false;
    }

    Layer* layer_buf = NULL;
    if ( bMoveUp )
    {
      layer_buf = unlocked_layers[0];
      for ( int i = 1; i < unlocked_layers.size(); i++ )
      {
        unlocked_layers[i-1] = unlocked_layers[i];
      }
      unlocked_layers[unlocked_layers.size()-1] = layer_buf;
    }
    else
    {
      layer_buf = unlocked_layers[unlocked_layers.size()-1];
      for ( int i = unlocked_layers.size()-1; i >= 1; i-- )
      {
        unlocked_layers[i] = unlocked_layers[i-1];
      }
      unlocked_layers[0] = layer_buf;
    }

    // put cycled unlocked layers back
    for ( int i = 0; i < m_layers.size(); i++ )
    {
      if ( m_layers[i]->IsLocked() )
      {
        if ( i < unlocked_layers.size() )
        {
          unlocked_layers.insert( unlocked_layers.begin() + i, m_layers[i] );
        }
        else
        {
          unlocked_layers.push_back( m_layers[i] );
        }
      }
    }
    m_layers = unlocked_layers;

    // restore visibility
    for ( int i = 0; i < m_layers.size(); i++ )
    {
      m_layers[i]->SetVisible( bVisibility[i] );
    }

    delete[] bVisibility;

    if ( nActive >= 0 )
    {
      SetActiveLayer( m_layers[nActive] );
    }

    emit LayerCycled( layer_buf );
    emit LayerMoved ( layer_buf );

    return true;
  }
  else
  {
    return false;
  }
}

bool LayerCollection::Contains( Layer* layer )
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i] == layer )
    {
      return true;
    }
  }
  return false;
}

void LayerCollection::Append2DProps( vtkRenderer* renderer, int nImagePlane )
{
  for ( int i = (int)m_layers.size()-1; i >= 0; i--)
  {
    m_layers[i]->Append2DProps( renderer, nImagePlane );
  }
}

void LayerCollection::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for ( int i = (int)m_layers.size()-1; i >= 0; i-- )
  {
    m_layers[i]->Append3DProps( renderer, bSliceVisibility );
  }
}

Layer* LayerCollection::GetFirstVisibleLayer()
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i]->IsVisible() )
    {
      return m_layers[i];
    }
  }
  return NULL;
}

int LayerCollection::GetNumberOfLayers()
{
  return (int)m_layers.size();
}

Layer* LayerCollection::GetLayer( int n )
{
  if ( n < (int)m_layers.size() )
  {
    return m_layers[n];
  }
  else
  {
    return NULL;
  }
}

void LayerCollection::SetActiveLayer( Layer* layer )
{
  if ( layer == NULL || this->Contains( layer ) )
  {
    m_layerActive = layer;
    emit ActiveLayerChanged( layer );
  }
}

Layer* LayerCollection::GetActiveLayer()
{
  return m_layerActive;
}

bool LayerCollection::SetSlicePosition( int nPlane, double dPos_in, bool bRoundToGrid )
{
  double dPos = dPos_in;

  if ( bRoundToGrid )
  {
    dPos = ((int)( ( dPos - m_dWorldOrigin[nPlane]) / m_dWorldVoxelSize[nPlane] ) ) * m_dWorldVoxelSize[nPlane]
           + m_dWorldOrigin[nPlane];
    if ( m_dSlicePosition[nPlane] <= m_dWorldOrigin[nPlane] + m_dWorldSize[nPlane] &&
         m_dSlicePosition[nPlane] >= m_dWorldOrigin[nPlane] &&
         ( dPos >  m_dWorldOrigin[nPlane] + m_dWorldSize[nPlane] || dPos < m_dWorldOrigin[nPlane] ) )
    {
      return false;
    }
  }

  if ( fabs( dPos - m_dSlicePosition[nPlane] ) < 1e-8 )
  {
    return false;
  }

  m_dSlicePosition[nPlane] = dPos;
  this->blockSignals( true );
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    m_layers[i]->SetSlicePosition( nPlane, dPos );
  }
  this->blockSignals( false );
  emit LayerActorUpdated();

  return true;
}

bool LayerCollection::OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid )
{
  return SetSlicePosition( nPlane, m_dSlicePosition[nPlane] + dPosDiff, bRoundToGrid );
}

bool LayerCollection::SetSlicePosition( double* slicePos )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = slicePos[i];
  }

  this->blockSignals( true );
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    m_layers[i]->SetSlicePosition( slicePos );
  }
  this->blockSignals( false );
  emit LayerActorUpdated();

  return true;
}

bool LayerCollection::SetSlicePosition( int nPlane, int nSliceNumber )
{
  return true;
}

double* LayerCollection::GetSlicePosition()
{
  return m_dSlicePosition;
}

void LayerCollection::GetSlicePosition( double* slicePos )
{
  for ( int i = 0; i < 3; i++ )
  {
    slicePos[i] = m_dSlicePosition[i];
  }
}

double* LayerCollection::GetCurrentRASPosition()
{
  return m_dCurrentRASPosition;
}

void LayerCollection::GetCurrentRASPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = m_dCurrentRASPosition[i];
  }
}

void LayerCollection::SetCurrentRASPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dCurrentRASPosition[i] = pos[i];
  }

  emit MouseRASPositionChanged();
}

double* LayerCollection::GetCursorRASPosition()
{
  return m_dCursorRASPosition;
}

void LayerCollection::GetCursorRASPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = m_dCursorRASPosition[i];
  }
}

void LayerCollection::SetCursorRASPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dCursorRASPosition[i] = pos[i];
  }

  emit CursorRASPositionChanged();
}

void LayerCollection::GetCurrentRASIndex( int* nIdx )
{
  for ( int i = 0; i < 3; i++ )
  {
    nIdx[i] = m_nCurrentRASIndex[i];
  }
}

QList<Layer*> LayerCollection::GetLayers()
{
  return m_layers;
}

QString LayerCollection::GetType()
{
  return m_strType;
}

double* LayerCollection::GetWorldOrigin()
{
  return m_dWorldOrigin;
}

void LayerCollection::SetWorldOrigin( double* dWorldOrigin )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldOrigin[i] = dWorldOrigin[i];
  }
}

double* LayerCollection::GetWorldSize()
{
  return m_dWorldSize;
}

void LayerCollection::SetWorldSize( double* dWorldSize )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldSize[i] = dWorldSize[i];
  }
}

double* LayerCollection::GetWorldVoxelSize()
{
  return m_dWorldVoxelSize;
}

void LayerCollection::SetWorldVoxelSize( double* dVoxelSize )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldVoxelSize[i] = dVoxelSize[i];
  }
}

void LayerCollection::GetWorldCenter( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = ( m_dWorldSize[i] + m_dWorldOrigin[i] ) / 2;
  }
}

Layer* LayerCollection::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < m_layers.size(); i++ )
  {
    if ( m_layers[i]->HasProp( prop ) )
    {
      return m_layers[i];
    }
  }
  return NULL;
}

void LayerCollection::LockCurrent( bool bLock )
{
  if ( m_layerActive )
  {
    m_layerActive->Lock( bLock );
  }
}

Layer* LayerCollection::GetLayer(const QString& type)
{
  for (int i = 0; i < m_layers.size(); i++)
  {
    if (m_layers[i]->IsTypeOf(type))
      return m_layers[i];
  }
  return NULL;
}

QList<Layer*> LayerCollection::GetLayers(const QString& type)
{
  QList<Layer*> layers;
  for (int i = 0; i < m_layers.size(); i++)
  {
    if (m_layers[i]->IsTypeOf(type))
      layers << m_layers[i];
  }
  return layers;
}
