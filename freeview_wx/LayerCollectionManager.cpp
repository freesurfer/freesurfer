/**
 * @file  LayerCollectionManager.cpp
 * @brief Manage the collections of layers, such as volumes, surfaces and way points, etc.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
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

#include "LayerCollectionManager.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include <math.h>

LayerCollectionManager::LayerCollectionManager() :
    Listener( "LayerCollectionManager" ),
    Broadcaster( "LayerCollectionManager" )
{
  LayerCollection* lc = new LayerCollection( "MRI" );
  lc->AddListener( this );
  m_layerCollections.push_back( lc );

  lc = new LayerCollection( "ROI" );
  lc->AddListener( this );
  m_layerCollections.push_back( lc );

  lc = new LayerCollection( "Surface" );
  lc->AddListener( this );
  m_layerCollections.push_back( lc );

  lc = new LayerCollection( "WayPoints" );
  lc->AddListener( this );
  m_layerCollections.push_back( lc );
}

LayerCollectionManager::~LayerCollectionManager()
{
  // safer way to delete all layer collections
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    m_layerCollections[i]->ClearAll( );
  }
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    delete m_layerCollections[i];
  }
  m_layerCollections.clear();
}

std::vector<Layer*> LayerCollectionManager::GetAllLayers()
{
  std::vector<Layer*> layers;

  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    std::vector<Layer*> sublayers = m_layerCollections[i]->GetLayers();
    for ( size_t j = 0; j < sublayers.size(); j++ )
      layers.push_back( sublayers[j] );
  }

  return layers;
}

void LayerCollectionManager::Append2DProps( vtkRenderer* renderer, int nImagePlane )
{
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    m_layerCollections[i]->Append2DProps( renderer, nImagePlane );
  }
}

void LayerCollectionManager::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    m_layerCollections[i]->Append3DProps( renderer, bSliceVisibility );
  }
}

void LayerCollectionManager::DoListenToMessage( std::string const iMsg, void* iData, void* sender )
{
// if ( iMsg == "LayerActorUpdated" )
  this->SendBroadcast( iMsg, iData, sender );
}

LayerCollection* LayerCollectionManager::GetLayerCollection( std::string strType )
{
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    if ( m_layerCollections[i]->GetType() == strType )
      return m_layerCollections[i];
  }
  return NULL;
}

bool LayerCollectionManager::SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid )
{
  bool bRet = false;
  this->BlockBroadcast( true );
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    if ( m_layerCollections[i]->SetSlicePosition( nPlane, dPos, bRoundToGrid ) )
      bRet = true;
  }
  this->BlockBroadcast( false );
  if ( bRet )
  {
    this->SendBroadcast( "LayerActorUpdated", this );
    this->SendBroadcast( "SlicePositionChanged", this );
  }

  return bRet;
}

bool LayerCollectionManager::SetSlicePosition( double* pos )
{
  bool bRet = false;
  this->BlockBroadcast( true );
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    if ( m_layerCollections[i]->SetSlicePosition( pos ) )
      bRet = true;
  }
  this->BlockBroadcast( false );
  if ( bRet )
  {
    this->SendBroadcast( "LayerActorUpdated", this );
    this->SendBroadcast( "SlicePositionChanged", this );
  }

  return bRet;
}

bool LayerCollectionManager::OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid  )
{
  bool bRet = false;
  this->BlockBroadcast( true );
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    if ( m_layerCollections[i]->OffsetSlicePosition( nPlane, dPosDiff, bRoundToGrid ) )
      bRet = true;
  }
  this->BlockBroadcast( false );
  if ( bRet )
  {
    this->SendBroadcast( "LayerActorUpdated", this );
    this->SendBroadcast( "SlicePositionChanged", this );
  }

  return bRet;
}

bool LayerCollectionManager::HasAnyLayer()
{
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    if ( !m_layerCollections[i]->IsEmpty() )
      return true;
  }

  return false;
}

bool LayerCollectionManager::HasLayer( std::string type )
{
  LayerCollection* lc = GetLayerCollection( type );

  return lc && !lc->IsEmpty();
}

void LayerCollectionManager::RefreshSlices()
{
  this->BlockBroadcast( true );
  for ( size_t i = 0; i < m_layerCollections.size(); i++ )
  {
    for ( int j = 0; j < m_layerCollections[i]->GetNumberOfLayers(); j++ )
    {
      for ( int n = 0; n < 3; n++ )
        m_layerCollections[i]->GetLayer( j )->OnSlicePositionChanged( n );
    }
  }

  this->BlockBroadcast( false );
  this->SendBroadcast( "LayerActorUpdated", this );
}
