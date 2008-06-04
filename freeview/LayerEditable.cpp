/**
 * @file  LayerEditable.cpp
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.4.2.1 $
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
 

#include <wx/wx.h>
#include "LayerEditable.h"
#include "vtkImageData.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "BrushProperty.h"
#include "MainWindow.h"
#include <stdlib.h>

LayerEditable::LayerEditable() : Layer()
{
	m_strTypeNames.push_back( "Editable" );
	
	m_fFillValue = 1;
	m_fBlankValue = 0;
	m_nBrushRadius = 1;
	m_bModified = false;
	m_bEditable = true;
	m_nMaxUndoSteps = 100;
	m_bufferClipboard.data = NULL;
	m_nActiveFrame = 0;
	m_propertyBrush = MainWindow::GetMainWindowPointer()->GetBrushProperty();
}

LayerEditable::~LayerEditable()
{
	for ( int i = 0; i < (int)m_bufferUndo.size(); i++ )
		delete[] m_bufferUndo[i].data;
	
	m_bufferUndo.clear();
	
	for ( int i = 0; i < (int)m_bufferRedo.size(); i++ )
		delete[] m_bufferRedo[i].data;
	
	if ( m_bufferClipboard.data )
		delete[] m_bufferClipboard.data;
	
	m_bufferUndo.clear();
}

bool LayerEditable::SetVoxelByIndex( int* n_in, int nPlane, bool bAdd )
{
	int* nDim = m_volumeRAS->GetDimensions();
/*	for ( int i = 0; i < 3; i++ )
	{
		if ( n_in[i] < 0 || n_in[i] >= nDim[i] )
			return false;
	}
*/	
//	float* ptr = ( float* )m_volumeRAS->GetScalarPointer( n );
//	if ( !ptr )
//		return false;
	
	int nBrushSize = m_propertyBrush->GetBrushSize();	
	int n[3], nsize[3] = { nBrushSize/2+1, nBrushSize/2+1, nBrushSize/2+1 };
	nsize[nPlane] = 1;
	int nActiveComp = GetActiveFrame();
	for ( int i = -nsize[0]+1; i < nsize[0]; i++ )
	{
		for ( int j = -nsize[1]+1; j < nsize[1]; j++ )
		{
			for ( int k = -nsize[2]+1; k < nsize[2]; k++ )
			{
				n[0] = n_in[0] + i;
				n[1] = n_in[1] + j;
				n[2] = n_in[2] + k;
				if ( n[0] >= 0 && n[0] < nDim[0] && 
					 n[1] >= 0 && n[1] < nDim[1] && 
					 n[2] >= 0 && n[2] < nDim[2] && 
					 MyUtils::GetDistance<int>( n, n_in ) <= nBrushSize/2.0 )
				{
				//	float fvalue = m_volumeRAS->GetScalarComponentAsFloat( n[0], n[1], n[2], 0 );
					if ( bAdd )
					{
						m_volumeRAS->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fFillValue );
					}
					else 
					{
						m_volumeRAS->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fBlankValue );
					}
				}
			}
		}
	}
	return true;
}

void LayerEditable::SetVoxelByRAS( double* ras, int nPlane, bool bAdd )
{
	int n[3];
	for ( int i = 0; i < 3; i++ )
		n[i] = ( int )( ( ras[i] - m_dWorldOrigin[i] ) / m_dWorldVoxelSize[i] + 0.5 );
	
	if ( SetVoxelByIndex( n, nPlane, bAdd ) )
	{
		SetModified();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else
	{
	//	PopUndo();		// roi not really changed, so pop the previously saved undo step
	}
}

void LayerEditable::SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd )
{
	int n1[3], n2[3];
	for ( int i = 0; i < 3; i++ )
	{
		n1[i] = ( int )( ( ras1[i] - m_dWorldOrigin[i] ) / m_dWorldVoxelSize[i] + 0.5 );
		n2[i] = ( int )( ( ras2[i] - m_dWorldOrigin[i] ) / m_dWorldVoxelSize[i] + 0.5 );
	}
	
	if ( SetVoxelByIndex( n1, n2, nPlane, bAdd ) )
	{
		SetModified();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else
	{
	//	PopUndo();
	}
}
	
bool LayerEditable::SetVoxelByIndex( int* n1, int* n2, int nPlane, bool bAdd )
{
	int nx = 1, ny = 2;
	if ( nPlane == 1 )
	{
		nx = 0;
		ny = 2;
	}
	else if (  nPlane == 2 )
	{
		nx = 0;
		ny = 1;
	}
	int x0 = n1[nx], y0 = n1[ny], x1 = n2[nx], y1 = n2[ny];
	
	int dx = x1 - x0;
	int dy = y1 - y0;
	double t = 0.5;
	int n[3];
	bool bChanged = SetVoxelByIndex( n1, nPlane, bAdd );	
	if ( abs( dx ) > abs( dy ) )
	{
		double m = (double) dy / (double) dx;
		t += y0;
		dx = ( dx < 0 ? -1 : 1 );
		m *= dx;
		while ( x0 != x1 )
		{
			x0 += dx;
			t += m;
			n[nx] = x0;
			n[ny] = (int) t;
			n[nPlane] = n1[nPlane];
			bChanged = SetVoxelByIndex( n, nPlane, bAdd );
		}
	}
	else
	{
		double m = (double) dx / (double) dy;
		t += x0;
		dy = ( dy < 0 ? -1 : 1 );
		m *= dy;
		while ( y0 != y1 )
		{
			y0 += dy;
			t += m;
			n[nx] = (int) t;
			n[ny] = y0;
			n[nPlane] = n1[nPlane];
			bChanged = SetVoxelByIndex( n, nPlane, bAdd );
		}
	}
	return true;
}

void LayerEditable::FloodFillByRAS( double* ras, int nPlane, bool bAdd )
{
	int n[3];
	for ( int i = 0; i < 3; i++ )
		n[i] = ( int )( ( ras[i] - m_dWorldOrigin[i] ) / m_dWorldVoxelSize[i] + 0.5 );
	
	if ( FloodFillByIndex( n, nPlane, bAdd ) )
	{
		SetModified();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
}

bool LayerEditable::FloodFillByIndex( int* n, int nPlane, bool bAdd )
{
	int* nDim = m_volumeRAS->GetDimensions();
	int nx = 0, ny = 0, x = 0, y = 0;
	switch ( nPlane )
	{
		case 0:
			nx = nDim[1];
			ny = nDim[2];
			x = n[1];
			y = n[2];
			break;
		case 1:
			nx = nDim[0];
			ny = nDim[2];
			x = n[0];
			y = n[2];
			break;
		case 2:
			nx = nDim[0];
			ny = nDim[1];
			x = n[0];
			y = n[1];
			break;
	}
	char** mask = MyUtils::AllocateMatrix<char>( ny, nx );
	int i, j;
	int nTolerance = m_propertyBrush->GetBrushTolerance();
	LayerEditable* ref_layer = m_propertyBrush->GetReferenceLayer();
	vtkImageData* ref = m_volumeRAS;
	int nActiveCompRef = 0;
	if ( ref_layer != NULL )
	{
		ref = ref_layer->GetRASVolume();
		nActiveCompRef = ref_layer->GetActiveFrame();
	}
	int nActiveComp = this->GetActiveFrame();
	float fVoxelValue = ref->GetScalarComponentAsFloat( n[0], n[1], n[2], 0 );
	switch ( nPlane )
	{
		case 0:
			if ( ref == m_volumeRAS )
			{
				for ( i = 0; i < nx; i++ )
				{
					for ( j = 0; j < ny; j++ )
					{
						mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveCompRef ) - fVoxelValue ) <= nTolerance
								? 1 : 0 );
					}
				}
			}
			else
			{
				for ( i = 0; i < nx; i++ )
				{
					for ( j = 0; j < ny; j++ )
					{
						mask[j][i] = ( m_volumeRAS->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveComp ) <= 0  &&
								fabs( ref->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveCompRef ) - fVoxelValue ) <= nTolerance
								? 1 : 0 );
					}
				}
			}
			break;
		case 1:
			if ( ref == m_volumeRAS )
			{
				for ( i = 0; i < nx; i++ )
				{
					for ( int j = 0; j < ny; j++ )
					{
						mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveCompRef ) - fVoxelValue ) <= nTolerance  
								? 1 : 0 );
					}
				}
			}
			else
			{			
				for ( i = 0; i < nx; i++ )
				{
					for ( int j = 0; j < ny; j++ )
					{
						mask[j][i] = ( m_volumeRAS->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveComp ) <= 0  &&
								fabs( ref->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveCompRef ) - fVoxelValue ) <= nTolerance  
								? 1 : 0 );
					}
				}
			}
			break;
		case 2:
			if ( ref == m_volumeRAS )
			{
				for ( i = 0; i < nx; i++ )
				{
					for ( j = 0; j < ny; j++ )
					{
						mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveCompRef ) - fVoxelValue ) <= nTolerance  
								? 1 : 0 );
					}
				}
			}
			else
			{
				for ( i = 0; i < nx; i++ )
				{
					for ( j = 0; j < ny; j++ )
					{
						mask[j][i] = ( m_volumeRAS->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveComp ) <= 0  &&
								fabs( ref->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveCompRef ) - fVoxelValue ) <= nTolerance  
								? 1 : 0 );
					}
				}
			}
			break;
	}
	
	MyUtils::FloodFill( mask, x, y, 0, 0, nx-1, ny-1, 2, 0 );	// unfill
	
	switch ( nPlane )
	{
		case 0:
			for ( i = 0; i < nx; i++ )
			{
				for ( j = 0; j < ny; j++ )
				{
					if ( mask[j][i] == 2 )
						m_volumeRAS->SetScalarComponentFromFloat( n[nPlane], i, j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
				}
			}
			break;
		case 1:
			for ( i = 0; i < nx; i++ )
			{
				for ( j = 0; j < ny; j++ )
				{
					if ( mask[j][i] == 2 )
						m_volumeRAS->SetScalarComponentFromFloat( i, n[nPlane], j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
				}
			}
			break;
		case 2:
			for ( i = 0; i < nx; i++ )
			{
				for ( j = 0; j < ny; j++ )
				{
					if ( mask[j][i] == 2 )
						m_volumeRAS->SetScalarComponentFromFloat( i, j, n[nPlane], nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
				}
			}
			break;
	}
	
	MyUtils::FreeMatrix( mask, ny );
	
	return true;
}

bool LayerEditable::HasUndo()
{
	return m_bufferUndo.size() > 0;
}

bool LayerEditable::HasRedo()
{
	return m_bufferRedo.size() > 0;
}

void LayerEditable::Undo()
{
	if ( m_bufferUndo.size() > 0 )
	{
		UndoRedoBufferItem item = m_bufferUndo[m_bufferUndo.size()-1];
		m_bufferUndo.pop_back();
		
		UndoRedoBufferItem item2;
		SaveBufferItem( item2, item.plane, item.slice );
		m_bufferRedo.push_back( item2 );
		
		LoadBufferItem( item );
		delete[] item.data;
		
		SetModified();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
}

void LayerEditable::Redo()
{
	if ( m_bufferRedo.size() > 0 )
	{
		UndoRedoBufferItem item = m_bufferRedo[m_bufferRedo.size()-1];
		m_bufferRedo.pop_back();
		
		UndoRedoBufferItem item2;
		SaveBufferItem( item2, item.plane, item.slice );
		m_bufferUndo.push_back( item2 );
		
		LoadBufferItem( item );
		delete[] item.data;
		
		SetModified();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
}

void LayerEditable::SaveForUndo( int nPlane )
{
	int nSlice = ( int )( ( m_dSlicePosition[nPlane] - m_dWorldOrigin[nPlane] ) / m_dWorldVoxelSize[nPlane] + 0.5 );
	
	if ( (int)m_bufferUndo.size() >= m_nMaxUndoSteps )
	{
		delete[] m_bufferUndo[0].data;
		m_bufferUndo.erase( m_bufferUndo.begin() );
	} 
	
	UndoRedoBufferItem item;
	SaveBufferItem( item, nPlane, nSlice );
	m_bufferUndo.push_back( item );
	
	// clear redo buffer
	for ( int i = 0; i < (int)m_bufferRedo.size(); i++ )
		delete[] m_bufferRedo[i].data;
	m_bufferRedo.clear();
}

bool LayerEditable::IsValidToPaste( int nPlane )
{
	return ( m_bufferClipboard.data != NULL && m_bufferClipboard.plane == nPlane );
}

void LayerEditable::Copy( int nPlane )
{
	int nSlice = ( int )( ( m_dSlicePosition[nPlane] - m_dWorldOrigin[nPlane] ) / m_dWorldVoxelSize[nPlane] + 0.5 );
	SaveBufferItem( m_bufferClipboard, nPlane, nSlice );
}

void LayerEditable::Paste( int nPlane )
{
	SaveForUndo( nPlane );
	
	int nSlice = ( int )( ( m_dSlicePosition[nPlane] - m_dWorldOrigin[nPlane] ) / m_dWorldVoxelSize[nPlane] + 0.5 );
	m_bufferClipboard.slice = nSlice;
	LoadBufferItem( m_bufferClipboard );
		
	SetModified();
	this->SendBroadcast( "LayerActorUpdated", this );
}

void LayerEditable::SaveBufferItem( UndoRedoBufferItem& item, int nPlane, int nSlice )
{
	int nDim[3], nStart[3] = { 0, 0, 0 };
	m_volumeRAS->GetDimensions( nDim );
	nDim[nPlane] = 1;
	nStart[nPlane] = nSlice;
	item.plane = nPlane;
	item.slice = nSlice;
	item.data = new char[nDim[0]*nDim[1]*nDim[2]*m_volumeRAS->GetScalarSize()];
	for ( int i = nStart[0]; i < nStart[0] + nDim[0]; i++ )
	{
		for ( int j = nStart[1]; j < nStart[1] + nDim[1]; j++ )
		{
			for ( int k = nStart[2]; k < nStart[2] + nDim[2]; k++ )
			{
				memcpy( item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_volumeRAS->GetScalarSize(),
						m_volumeRAS->GetScalarPointer( i, j, k ), 
						m_volumeRAS->GetScalarSize() );
			}
		}
	}
}
		
void LayerEditable::LoadBufferItem( UndoRedoBufferItem& item )
{
	int nDim[3], nStart[3] = { 0, 0, 0 };
	m_volumeRAS->GetDimensions( nDim );
	nDim[item.plane] = 1;
	nStart[item.plane] = item.slice;
	for ( int i = nStart[0]; i < nStart[0] + nDim[0]; i++ )
	{
		for ( int j = nStart[1]; j < nStart[1] + nDim[1]; j++ )
		{
			for ( int k = nStart[2]; k < nStart[2] + nDim[2]; k++ )
			{
				memcpy( m_volumeRAS->GetScalarPointer( i, j, k ),
						item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_volumeRAS->GetScalarSize(),
						m_volumeRAS->GetScalarSize() );						
			}
		}
	}
}

float LayerEditable::GetFillValue()
{
	return m_fFillValue;
}

void LayerEditable::SetFillValue( float fFill )
{
	m_fFillValue = fFill;
}
		
float LayerEditable::GetBlankValue()
{
	return m_fBlankValue;
}

void LayerEditable::SetBlankValue( float fBlank )
{
	m_fBlankValue = fBlank;
}

int LayerEditable::GetBrushRadius()
{
	return m_nBrushRadius;
}

void LayerEditable::SetBrushRadius( int nRadius )
{
	m_nBrushRadius = nRadius;
}
		
void LayerEditable::SetModified()
{
	m_bModified = true;
	this->SendBroadcast( "LayerModified", this );
}
