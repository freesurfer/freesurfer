/**
 * @file  LayerCollection.cpp
 * @brief A collection of layers.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/06/11 21:30:18 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
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
 
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <math.h>

LayerCollection::LayerCollection( std::string strType) : 
	Listener( "LayerCollection" ), 
	Broadcaster( "LayerCollection" ),
	m_layerActive( 0 ),
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
	for ( unsigned int i = 0; i < m_layers.size(); i++ )
		delete m_layers[i];

	m_layers.clear();
}
 
bool LayerCollection::IsEmpty()
{
	return m_layers.size() == 0;
}

int LayerCollection::GetLayerIndex( Layer* layer )
{
	for ( unsigned int i = 0; i < m_layers.size(); i++ )
	{
		if ( m_layers[i] == layer )
			return i;
	}
	return -1;
}
 
bool LayerCollection::AddLayer( Layer* layer, bool initializeCoordinate )
{
	if ( !layer->IsTypeOf( m_strType ) )
	{
		std::cerr << "Can not add layer type of " << layer->GetEndType().c_str() 
				<< " to layer collection type of " <<  m_strType.c_str() << endl;
		return false;
	}
	
	for ( unsigned int i = 0; i < m_layers.size(); i++ )
	{
		if ( m_layers[i] == layer )
			return false;
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
	layer->AddListener( this );
	
	this->SetActiveLayer( layer );
	
	this->SendBroadcast( "LayerAdded", layer );
	
	return true;
}

bool LayerCollection::RemoveLayer( Layer* layer, bool deleteObject )
{
	for ( unsigned int i = 0; i < m_layers.size(); i++ )
	{
		if ( m_layers[i] == layer )
		{
			if (deleteObject)
				delete m_layers[i];
			
			m_layers.erase( m_layers.begin() + i );
			
			if ( m_layers.size() == 0 )
				SetActiveLayer( NULL );
			else 
			{
				if ( i == m_layers.size() )
					SetActiveLayer( m_layers[m_layers.size()-1] );
				else
					SetActiveLayer( m_layers[i] );
			} 
			
			this->SendBroadcast( "LayerRemoved", layer );
			
			return true;
		}
	}
	
	return false;
}

bool LayerCollection::MoveLayerUp( Layer* layer )
{
	for ( unsigned int i = 1; i < m_layers.size(); i++)
	{
		if ( m_layers[i] == layer )
		{
			Layer* temp = m_layers[i-1];
			m_layers[i-1] = layer;
			m_layers[i] = temp;
			
			this->SendBroadcast( "LayerMoved", layer );
			
			return true;
		}
	}
	return false;
}

bool LayerCollection::MoveLayerDown( Layer* layer )
{
	for ( unsigned int i = 0; i < m_layers.size()-1; i++)
	{
		if ( m_layers[i] == layer )
		{
			Layer* temp = m_layers[i+1];
			m_layers[i+1] = layer;
			m_layers[i] = temp;
			
			this->SendBroadcast( "LayerMoved", layer );
			
			return true;
		}
	}
	return false;
}

bool LayerCollection::CycleLayer()
{
	if ( (int)m_layers.size() > 1 )
	{
		Layer* layer0 = m_layers[0];
		for ( unsigned int i = 1; i < m_layers.size(); i++ )
		{
			m_layers[i-1] = m_layers[i];
		}
		m_layers[m_layers.size()-1] = layer0;
		
		this->SendBroadcast( "LayerCycled", layer0 );
		this->SendBroadcast( "LayerMoved", layer0 );
		
		return true;
	}
	else
		return false;
}

bool LayerCollection::Contains( Layer* layer )
{
	for ( unsigned int i = 0; i < m_layers.size(); i++)
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

void LayerCollection::Append3DProps( vtkRenderer* renderer )
{
	for ( int i = (int)m_layers.size()-1; i >= 0; i-- )
	{
		m_layers[i]->Append3DProps( renderer );
	}
}

int LayerCollection::GetNumberOfLayers()
{
	return (int)m_layers.size();
}

Layer* LayerCollection::GetLayer( int n )
{
	return m_layers[n];
}

void LayerCollection::DoListenToMessage( std::string const iMsg, void* iData )
{
//	if ( iMsg == "LayerActorUpdated" )
		this->SendBroadcast( iMsg, iData );
}

void LayerCollection::SetActiveLayer( Layer* layer )
{
	if ( layer == NULL || this->Contains( layer ) )
		m_layerActive = layer;
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
	}
	
	if ( dPos > m_dWorldOrigin[nPlane] + m_dWorldSize[nPlane] || 
			dPos < m_dWorldOrigin[nPlane] ||
			fabs( dPos - m_dSlicePosition[nPlane] ) < 1e-8 )
		return false;
	
	m_dSlicePosition[nPlane] = dPos;
	this->BlockBroadcast( true );
	for ( int i = 0; i < (int)m_layers.size(); i++ )
	{
		m_layers[i]->SetSlicePosition( nPlane, dPos );
	}
	this->BlockBroadcast( false );
	this->SendBroadcast( "LayerActorUpdated", this );
	
	return true;
}

bool LayerCollection::OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid )
{
	return SetSlicePosition( nPlane, m_dSlicePosition[nPlane] + dPosDiff, bRoundToGrid );
}
		
bool LayerCollection::SetSlicePosition( double* slicePos )
{
	for ( int i = 0; i < 3; i++ )
		m_dSlicePosition[i] = slicePos[i];
	
	this->BlockBroadcast( true );
	for ( int i = 0; i < (int)m_layers.size(); i++ )
	{
		m_layers[i]->SetSlicePosition( slicePos );
	}
	this->BlockBroadcast( false );
	this->SendBroadcast( "LayerActorUpdated", this );
	
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
		slicePos[i] = m_dSlicePosition[i];
}

double* LayerCollection::GetCurrentRASPosition()
{
	return m_dCurrentRASPosition;
}

void LayerCollection::GetCurrentRASPosition( double* pos )
{
	for ( int i = 0; i < 3; i++ )
		pos[i] = m_dCurrentRASPosition[i];
}

void LayerCollection::SetCurrentRASPosition( double* pos )
{
	for ( int i = 0; i < 3; i++ )
		m_dCurrentRASPosition[i] = pos[i];
	
	this->SendBroadcast( "MouseRASPositionChanged", this );
}
		

double* LayerCollection::GetCursorRASPosition()
{
	return m_dCursorRASPosition;
}

void LayerCollection::GetCursorRASPosition( double* pos )
{
	for ( int i = 0; i < 3; i++ )
		pos[i] = m_dCursorRASPosition[i];
}

void LayerCollection::SetCursorRASPosition( double* pos )
{
	for ( int i = 0; i < 3; i++ )
		m_dCursorRASPosition[i] = pos[i];
	
	this->SendBroadcast( "CursorRASPositionChanged", this );
}
			
void LayerCollection::GetCurrentRASIndex( int* nIdx )
{
	for ( int i = 0; i < 3; i++ )
		nIdx[i] = m_nCurrentRASIndex[i];
}

std::vector<Layer*> LayerCollection::GetLayers()
{
	return m_layers;
}

std::string LayerCollection::GetType()
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
		m_dWorldOrigin[i] = dWorldOrigin[i];
}
		
double* LayerCollection::GetWorldSize()
{
	return m_dWorldSize;
}	

void LayerCollection::SetWorldSize( double* dWorldSize )
{
	for ( int i = 0; i < 3; i++ )
		m_dWorldSize[i] = dWorldSize[i];
}
		
double* LayerCollection::GetWorldVoxelSize()
{
	return m_dWorldVoxelSize;
}

void LayerCollection::SetWorldVoxelSize( double* dVoxelSize )
{
	for ( int i = 0; i < 3; i++ )
		m_dWorldVoxelSize[i] = dVoxelSize[i];
}

bool LayerCollection::HasProp( vtkProp* prop )
{
	for ( int i = 0; i < (int)m_layers.size(); i++ )
	{
		if ( m_layers[i]->HasProp( prop ) )
			return true;
	}
	return false;
}
