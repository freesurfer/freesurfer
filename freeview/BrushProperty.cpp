/**
 * @file  BrushProperty.cpp
 * @brief class to hold brush properties
 *
 * Simple mix-in class for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/07 22:01:54 $
 *    $Revision: 1.5 $
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


#include "BrushProperty.h"
#include "LayerVolumeBase.h"
#include <wx/config.h>

using namespace std;

BrushProperty::BrushProperty () :
		m_nBrushSize( 1 ), 
		m_nBrushTolerance( 0 ), 
		m_bEnableDrawRange( false ),
		m_bEnableExcludeRange( false ),
		m_bDrawConnectedOnly( false ),
		m_layerRef( NULL )
{
	m_dDrawRange[0] = 0;
	m_dDrawRange[1] = 1000000;
	m_dExcludeRange[0] = 0;
	m_dExcludeRange[1] = 0;
	wxConfigBase* config = wxConfigBase::Get();
	if ( config )
	{
		config->Read( _T("/BrushProperty/Size"), &m_nBrushSize, 1L );
	//	config->Read( _T("/BrushProperty/Tolerance"), &m_nBrushTolerance, 0L );
		m_nBrushTolerance = 0;		
		config->Read( _T("/BrushProperty/EnableDrawRange"), &m_bEnableDrawRange, false );
		config->Read( _T("/BrushProperty/EnableExcludeRange"), &m_bEnableExcludeRange, false );
		config->Read( _T("/BrushProperty/DrawConnected"), &m_bDrawConnectedOnly, false );
		config->Read( _T("/BrushProperty/DrawRangeLow"), m_dDrawRange, 0 );
		config->Read( _T("/BrushProperty/DrawRangeHigh"), m_dDrawRange+1, 1000000 );
		config->Read( _T("/BrushProperty/ExcludeRangeLow"), m_dExcludeRange, 0 );
		config->Read( _T("/BrushProperty/ExcludeRangeHigh"), m_dExcludeRange+1, 0 );
	}
}
	
BrushProperty::~BrushProperty()
{
	wxConfigBase* config = wxConfigBase::Get();
	if ( config )
	{
		config->Write( _T("/BrushProperty/Size"), m_nBrushSize );
		config->Write( _T("/BrushProperty/Tolerance"), m_nBrushTolerance );
		config->Write( _T("/BrushProperty/EnableDrawRange"), m_bEnableDrawRange );
		config->Write( _T("/BrushProperty/EnableExcludeRange"), m_bEnableExcludeRange );
		config->Write( _T("/BrushProperty/DrawConnected"), m_bDrawConnectedOnly );
		config->Write( _T("/BrushProperty/DrawRangeLow"), m_dDrawRange[0] );
		config->Write( _T("/BrushProperty/DrawRangeHigh"), m_dDrawRange[1] );
		config->Write( _T("/BrushProperty/ExcludeRangeLow"), m_dExcludeRange[0] );
		config->Write( _T("/BrushProperty/ExcludeRangeHigh"), m_dExcludeRange[1] );
	}	
}
	
int BrushProperty::GetBrushSize()
{
	return m_nBrushSize;
}

void BrushProperty::SetBrushSize( int nSize )
{
	m_nBrushSize = nSize;
}
	
int BrushProperty::GetBrushTolerance()
{
	return m_nBrushTolerance;
}

void BrushProperty::SetBrushTolerance( int nTolerance )
{
	m_nBrushTolerance = nTolerance;
}

LayerVolumeBase* BrushProperty::GetReferenceLayer()
{
	return m_layerRef;
}
	
void BrushProperty::SetReferenceLayer( LayerVolumeBase* layer )
{
	m_layerRef = layer;
}

double* BrushProperty::GetDrawRange()
{
	return m_dDrawRange;
}
	
void BrushProperty::SetDrawRange( double* range )
{
	SetDrawRange( range[0], range[1] );
}
	
void BrushProperty::SetDrawRange( double low, double high )
{
	m_dDrawRange[0] = low;
	m_dDrawRange[1] = high;
}
	
	
bool BrushProperty::GetDrawRangeEnabled()
{
	return m_bEnableDrawRange;
}
	
void BrushProperty::SetDrawRangeEnabled( bool bEnable )
{
	m_bEnableDrawRange = bEnable;
} 
	 
double* BrushProperty::GetExcludeRange()
{
	return m_dExcludeRange;
}

void BrushProperty::SetExcludeRange( double* range )
{
	SetExcludeRange( range[0], range[1] );
}

void BrushProperty::SetExcludeRange( double low, double high )
{
	m_dExcludeRange[0] = low;
	m_dExcludeRange[1] = high;
}
	
bool BrushProperty::GetExcludeRangeEnabled()
{
	return m_bEnableExcludeRange;
}

void BrushProperty::SetExcludeRangeEnabled( bool bEnable )
{
	m_bEnableExcludeRange = bEnable;
} 

bool BrushProperty::GetDrawConnectedOnly()
{
	return m_bDrawConnectedOnly;
}

void BrushProperty::SetDrawConnectedOnly( bool bEnable )
{
	m_bDrawConnectedOnly = bEnable;
}
