/**
 * @file  LayerCollection.h
 * @brief LayerCollection data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.2.2.1 $
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
 
#ifndef LayerCollection_h
#define LayerCollection_h

#include <string>
#include <vector>
#include "Listener.h"
#include "Broadcaster.h"

class Layer;
class LayerMRI;
class vtkRenderer;

class LayerCollection : public Listener, public Broadcaster
{
	public:
		LayerCollection( std::string type );
		virtual ~LayerCollection();	
		
		int GetNumberOfLayers();		
		Layer* GetLayer( int i );
		
		bool AddLayer( Layer* layer, bool initializeCoordinate = false );
		bool RemoveLayer( Layer* layer, bool deleteObject = true );
		bool MoveLayerUp( Layer* layer );
		bool MoveLayerDown( Layer* layer );
    
		void Append2DProps( vtkRenderer* renderer, int nImagePlane );
		void Append3DProps( vtkRenderer* renderer );
		
		bool Contains( Layer* layer );
		bool IsEmpty();
		
		void SetActiveLayer( Layer* layer );
		Layer* GetActiveLayer();
		
		double* GetSlicePosition();
		void GetSlicePosition( double* slicePos );
		bool SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid = true );
		bool OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid = true );
		bool SetSlicePosition( int nPlane, int nSliceNumber );
		bool SetSlicePosition( double* slicePos );
		
		double* GetCurrentRASPosition();
		void GetCurrentRASPosition( double* pos );
		void SetCurrentRASPosition( double* pos );
		
		double* GetCursorRASPosition();
		void GetCursorRASPosition( double* pos );
		void SetCursorRASPosition( double* pos );
		
		void GetCurrentRASIndex( int* nIdx );

		double* GetWorldOrigin();
		void SetWorldOrigin( double* dWorldOrigin );
		
		double* GetWorldSize();		
		void SetWorldSize( double* dWorldSize );
		
		double* GetWorldVoxelSize();
		void SetWorldVoxelSize( double* dVoxelSize );	
		
		std::vector<Layer*> GetLayers();
		
		virtual void DoListenToMessage( std::string const iMsg, void* iData );
		
		std::string GetType();
				
	protected:
		std::vector<Layer*>	m_layers;
		
		double 	m_dSlicePosition[3];
		double	m_dWorldOrigin[3];
		double	m_dWorldSize[3];		
		double	m_dWorldVoxelSize[3];
		
		double	m_dCurrentRASPosition[3];
		int		m_nCurrentRASIndex[3];
		
		double	m_dCursorRASPosition[3];
		
		Layer*	m_layerActive;
		std::string		m_strType;
};

#endif 


