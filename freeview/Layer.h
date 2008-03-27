/**
 * @file  Layer.h
 * @brief Layer data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:15 $
 *    $Revision: 1.1 $
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
 
#ifndef Layer_h
#define Layer_h

#include "Listener.h"
#include "Broadcaster.h"
#include <string>
#include <vector>

class vtkRenderer;

class Layer : public Listener, public Broadcaster
{
	public:
		Layer();
		virtual ~Layer();	
		
		const char*	GetName()
			{ return m_strName.c_str(); }
		
		void SetName( const char* name )
			{ m_strName = name; }
    	
		virtual void Append2DProps( vtkRenderer* renderer, int nPlane ) = 0;
		virtual void Append3DProps( vtkRenderer* renderer ) = 0;
		
		virtual void SetVisible( bool bVisible = true ) = 0;
		virtual bool IsVisible() = 0;
		
		double* GetWorldOrigin();
		void GetWorldOrigin( double* origin );
		void SetWorldOrigin( double* origin );
		
		double* GetWorldVoxelSize();
		void GetWorldVoxelSize( double* voxelsize );
		void SetWorldVoxelSize( double* voxelsize );
		
		double* GetWorldSize();
		void GetWorldSize( double* size );
		void SetWorldSize( double* size );
		
		double* GetSlicePosition();
		void GetSlicePosition( double* slicePos );
		void SetSlicePosition( double* slicePos );
		void SetSlicePosition( int nPlane, double slicePos );
		
		virtual void OnSlicePositionChanged( int nPlane ) = 0;
		
		bool IsTypeOf( std::string tname );
		
		std::string GetErrorString()
			{ return m_strError; }
		
		void SetErrorString( std::string msg )
			{ m_strError = msg; }
		
		std::string GetEndType();
		
	protected:
		std::string		m_strName;
		double			m_dSlicePosition[3];
		double			m_dWorldOrigin[3];
		double			m_dWorldVoxelSize[3];
		double 			m_dWorldSize[3];
		std::string		m_strError;
		std::vector<std::string>	m_strTypeNames;	
};

#endif 


