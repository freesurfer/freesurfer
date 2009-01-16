/**
 * @file  LayerVolumeBase.cpp
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/01/16 22:13:07 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2009,
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
#include "LayerVolumeBase.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "BrushProperty.h"
#include "MainWindow.h"
#include "LivewireTool.h"
#include <stdlib.h>

LayerVolumeBase::LayerVolumeBase() : LayerEditable()
{
	m_strTypeNames.push_back( "VolumeBase" );
	
	m_fFillValue = 1;
	m_fBlankValue = 0;
	m_nBrushRadius = 1;
	m_bufferClipboard.data = NULL;
	m_nActiveFrame = 0;
	m_propertyBrush = MainWindow::GetMainWindowPointer()->GetBrushProperty();
	m_livewire = new LivewireTool();
	m_imageData = NULL;
	m_imageDataRef = NULL;
}

LayerVolumeBase::~LayerVolumeBase()
{
	for ( size_t i = 0; i < m_bufferUndo.size(); i++ )
		delete[] m_bufferUndo[i].data;
	
	m_bufferUndo.clear();
	
	for ( size_t i = 0; i < m_bufferRedo.size(); i++ )
		delete[] m_bufferRedo[i].data;
	
	if ( m_bufferClipboard.data )
		delete[] m_bufferClipboard.data;
	
	m_bufferUndo.clear();
	
	delete m_livewire;
}

bool LayerVolumeBase::SetVoxelByIndex( int* n_in, int nPlane, bool bAdd )
{
	int* nDim = m_imageData->GetDimensions();
/*	for ( int i = 0; i < 3; i++ )
	{
		if ( n_in[i] < 0 || n_in[i] >= nDim[i] )
			return false;
	}
*/	
//	float* ptr = ( float* )m_imageData->GetScalarPointer( n );
//	if ( !ptr )
//		return false;
	
	int nBrushSize = m_propertyBrush->GetBrushSize();	
	int n[3], nsize[3] = { nBrushSize/2+1, nBrushSize/2+1, nBrushSize/2+1 };
	nsize[nPlane] = 1;
	int nActiveComp = GetActiveFrame();
	double* draw_range = m_propertyBrush->GetDrawRange();
	double* exclude_range = m_propertyBrush->GetExcludeRange();
	LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
	vtkImageData* ref = m_imageData;
	int nActiveCompRef = 0;
	if ( ref_layer != NULL )
	{
		ref = ref_layer->GetImageData();
		nActiveCompRef = ref_layer->GetActiveFrame();
	}
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
					if ( bAdd ) 
					{
						double fvalue = ref->GetScalarComponentAsDouble( n[0], n[1], n[2], nActiveCompRef );
						if ( ( m_propertyBrush->GetDrawRangeEnabled() && 
							 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
							 ( m_propertyBrush->GetExcludeRangeEnabled() && 
							 ( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) || 
							 ( m_propertyBrush->GetDrawConnectedOnly() && 
							 ( !GetConnectedToOld( m_imageData, nActiveComp, n, nPlane ) ) ) )
							;
						else	 
							m_imageData->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fFillValue );
					}
					else 
					{
						m_imageData->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fBlankValue );
					}
				}
			}
		}
	}
	return true;
}

bool LayerVolumeBase::GetConnectedToOld( vtkImageData* img, int nFrame, int* n_in, int nPlane )
{
	int* nDim = img->GetDimensions();
	int nBounds[6] = { -1, 1, -1, 1, -1, 1 };
	nBounds[nPlane*2] = nBounds[nPlane*2+1] = 0;
	int n[3];
	for ( int i = nBounds[0]; i <= nBounds[1]; i++ )
	{
		for ( int j = nBounds[2]; j <= nBounds[3]; j++ )
		{
			for ( int k = nBounds[4]; k <= nBounds[5]; k++ )
			{
				n[0] = n_in[0] + i;
				n[1] = n_in[1] + j;
				n[2] = n_in[2] + k;
				if ( abs( i + j + k ) == 1 &&
					 n[0] >= 0 && n[0] < nDim[0] && 
					 n[1] >= 0 && n[1] < nDim[1] && 
					 n[2] >= 0 && n[2] < nDim[2] )
				{
					double fvalue = img->GetScalarComponentAsDouble( n[0], n[1], n[2], nFrame );
					if ( fvalue > 0 )
						return true;
				}
			}
		}
	}
	return false;
}

void LayerVolumeBase::SetVoxelByRAS( double* ras, int nPlane, bool bAdd )
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

void LayerVolumeBase::SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd )
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
	
bool LayerVolumeBase::SetVoxelByIndex( int* n1, int* n2, int nPlane, bool bAdd )
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

void LayerVolumeBase::FloodFillByRAS( double* ras, int nPlane, bool bAdd )
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

bool LayerVolumeBase::FloodFillByIndex( int* n, int nPlane, bool bAdd )
{
	int* nDim = m_imageData->GetDimensions();
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
	double* draw_range = m_propertyBrush->GetDrawRange();
	double* exclude_range = m_propertyBrush->GetExcludeRange();
	LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
	vtkImageData* ref = m_imageData;
	int nActiveCompRef = 0;
	if ( ref_layer != NULL )
	{
		ref = ref_layer->GetImageData();
		nActiveCompRef = ref_layer->GetActiveFrame();
	}
	int nActiveComp = this->GetActiveFrame();
	float fVoxelValue = ref->GetScalarComponentAsFloat( n[0], n[1], n[2], 0 );
	switch ( nPlane )
	{
		case 0:
			if ( ref == m_imageData )
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
						mask[j][i] = ( m_imageData->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveComp ) <= 0  &&
								fabs( ref->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveCompRef ) - fVoxelValue ) <= nTolerance
								? 1 : 0 );
					}
				}
			}
			break;
		case 1:
			if ( ref == m_imageData )
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
						mask[j][i] = ( m_imageData->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveComp ) <= 0  &&
								fabs( ref->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveCompRef ) - fVoxelValue ) <= nTolerance  
								? 1 : 0 );
					}
				}
			}
			break;
		case 2:
			if ( ref == m_imageData )
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
						mask[j][i] = ( m_imageData->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveComp ) <= 0  &&
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
					{
						double fvalue = ref->GetScalarComponentAsDouble( n[nPlane], i, j, nActiveCompRef );
						if ( ( m_propertyBrush->GetDrawRangeEnabled() && 
							( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
							( m_propertyBrush->GetExcludeRangeEnabled() && 
							( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) )
							;
						else
						{
							m_imageData->SetScalarComponentFromFloat( n[nPlane], i, j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
						}
					}
				}
			}
			break;
		case 1:
			for ( i = 0; i < nx; i++ )
			{
				for ( j = 0; j < ny; j++ )
				{
					if ( mask[j][i] == 2 )
					{
						double fvalue = ref->GetScalarComponentAsDouble( i, n[nPlane], j, nActiveCompRef );
						if ( ( m_propertyBrush->GetDrawRangeEnabled() && 
								 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
								 ( m_propertyBrush->GetExcludeRangeEnabled() && 
								 ( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) )
							;
						else
						{
							m_imageData->SetScalarComponentFromFloat( i, n[nPlane], j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
						}
					}
				}
			}
			break;
		case 2:
			for ( i = 0; i < nx; i++ )
			{
				for ( j = 0; j < ny; j++ )
				{
					if ( mask[j][i] == 2 )
					{
						double fvalue = ref->GetScalarComponentAsDouble( i, j, n[nPlane], nActiveCompRef );
						if ( ( m_propertyBrush->GetDrawRangeEnabled() && 
							( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
							( m_propertyBrush->GetExcludeRangeEnabled() && 
							( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) )
							;
						else
						{
							m_imageData->SetScalarComponentFromFloat( i, j, n[nPlane], nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
						}
					}
				}
			}
			break;
	}
	
	MyUtils::FreeMatrix( mask, ny );
	
	return true;
}

void LayerVolumeBase::SetLiveWireByRAS( double* pt1, double* pt2, int nPlane )
{
	int n1[3], n2[3];
	double* orig = m_imageData->GetOrigin();
	double* vxlsize = m_imageData->GetSpacing();
	for ( int i = 0; i < 3; i++ )
	{
		n1[i] = ( int )( ( pt1[i] - orig[i] ) / vxlsize[i] );
		n2[i] = ( int )( ( pt2[i] - orig[i] ) / vxlsize[i] );
	}
	
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkImageData* image = m_imageData;
	LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
	if ( ref_layer != NULL )
		image = ref_layer->GetImageData();
//	else if ( m_imageDataRef.GetPointer() != NULL )
	else if ( m_imageDataRef != NULL )
		image = m_imageDataRef;
	
	m_livewire->GetLivewirePoints( image, nPlane, n1[nPlane], pt1, pt2, pts );
	int n[3];
	for ( int i = 0; i < pts->GetNumberOfPoints(); i++ )
	{
		double* p = pts->GetPoint( i );
		n[0] = (int)( ( p[0] - orig[0] ) / vxlsize[0] + 0.5 );
		n[1] = (int)( ( p[1] - orig[1] ) / vxlsize[1] + 0.5);
		n[2] = (int)( ( p[2] - orig[2] ) / vxlsize[2] + 0.5 );
			
		SetVoxelByIndex( n, nPlane, true );
	}
	
	SetModified();
	this->SendBroadcast( "LayerActorUpdated", this );
}

std::vector<double> LayerVolumeBase::GetLiveWirePointsByRAS( double* pt1, double* pt2, int nPlane )
{
	vtkImageData* image = m_imageData;
	LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
	if ( ref_layer != NULL )
		image = ref_layer->GetImageData();
//	else if ( m_imageDataRef.GetPointer() != NULL )
	else if ( m_imageDataRef != NULL )
		image = m_imageDataRef;
	
	int n1[3], n2[3];
	double* orig = image->GetOrigin();
	double* vxlsize = image->GetSpacing();
	for ( int i = 0; i < 3; i++ )
	{
		n1[i] = ( int )( ( pt1[i] - orig[i] ) / vxlsize[i] );
		n2[i] = ( int )( ( pt2[i] - orig[i] ) / vxlsize[i] );
	}
	
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	
	m_livewire->GetLivewirePoints( image, nPlane, n1[nPlane], pt1, pt2, pts );	
	std::vector<double> pts_out;
	for ( int i = 1; i < pts->GetNumberOfPoints()-1; i++ )
	{
		double* p = pts->GetPoint( i );	
		pts_out.push_back( p[0] ); //*vxlsize[0] + orig[0] );
		pts_out.push_back( p[1] ); //*vxlsize[1] + orig[1] );
		pts_out.push_back( p[2] ); //*vxlsize[2] + orig[2] );
	}
	return pts_out;
}
	
bool LayerVolumeBase::SetLiveWireByIndex( int* n1, int* n2, int nPlane )
{
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
//	MyUtils::GetLivewirePoints( m_imageData, nPlane, n1[nPlane], n1, n2, pts );
	return true;
}


bool LayerVolumeBase::HasUndo()
{
	return m_bufferUndo.size() > 0;
}

bool LayerVolumeBase::HasRedo()
{
	return m_bufferRedo.size() > 0;
}

void LayerVolumeBase::Undo()
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
		this->SendBroadcast( "LayerEdited", this );
	}
}

void LayerVolumeBase::Redo()
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
		this->SendBroadcast( "LayerEdited", this );
	}
}

void LayerVolumeBase::SaveForUndo( int nPlane )
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
	for ( size_t i = 0; i < m_bufferRedo.size(); i++ )
		delete[] m_bufferRedo[i].data;
	m_bufferRedo.clear();
}

bool LayerVolumeBase::IsValidToPaste( int nPlane )
{
	return ( m_bufferClipboard.data != NULL && m_bufferClipboard.plane == nPlane );
}

void LayerVolumeBase::Copy( int nPlane )
{
	int nSlice = ( int )( ( m_dSlicePosition[nPlane] - m_dWorldOrigin[nPlane] ) / m_dWorldVoxelSize[nPlane] + 0.5 );
	SaveBufferItem( m_bufferClipboard, nPlane, nSlice );
}

void LayerVolumeBase::Paste( int nPlane )
{
	SaveForUndo( nPlane );
	
	int nSlice = ( int )( ( m_dSlicePosition[nPlane] - m_dWorldOrigin[nPlane] ) / m_dWorldVoxelSize[nPlane] + 0.5 );
	m_bufferClipboard.slice = nSlice;
	LoadBufferItem( m_bufferClipboard );
		
	SetModified();
	this->SendBroadcast( "LayerActorUpdated", this );
}

void LayerVolumeBase::SaveBufferItem( UndoRedoBufferItem& item, int nPlane, int nSlice )
{
	int nDim[3], nStart[3] = { 0, 0, 0 };
	m_imageData->GetDimensions( nDim );
	nDim[nPlane] = 1;
	nStart[nPlane] = nSlice;
	item.plane = nPlane;
	item.slice = nSlice;
	item.data = new char[nDim[0]*nDim[1]*nDim[2]*m_imageData->GetScalarSize()];
	for ( int i = nStart[0]; i < nStart[0] + nDim[0]; i++ )
	{
		for ( int j = nStart[1]; j < nStart[1] + nDim[1]; j++ )
		{
			for ( int k = nStart[2]; k < nStart[2] + nDim[2]; k++ )
			{
				memcpy( item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_imageData->GetScalarSize(),
						m_imageData->GetScalarPointer( i, j, k ), 
						m_imageData->GetScalarSize() );
			}
		}
	}
}
		
void LayerVolumeBase::LoadBufferItem( UndoRedoBufferItem& item )
{
	int nDim[3], nStart[3] = { 0, 0, 0 };
	m_imageData->GetDimensions( nDim );
	nDim[item.plane] = 1;
	nStart[item.plane] = item.slice;
	for ( int i = nStart[0]; i < nStart[0] + nDim[0]; i++ )
	{
		for ( int j = nStart[1]; j < nStart[1] + nDim[1]; j++ )
		{
			for ( int k = nStart[2]; k < nStart[2] + nDim[2]; k++ )
			{
				memcpy( m_imageData->GetScalarPointer( i, j, k ),
						item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_imageData->GetScalarSize(),
						m_imageData->GetScalarSize() );						
			}
		}
	}
}

float LayerVolumeBase::GetFillValue()
{
	return m_fFillValue;
}

void LayerVolumeBase::SetFillValue( float fFill )
{
	m_fFillValue = fFill;
}
		
float LayerVolumeBase::GetBlankValue()
{
	return m_fBlankValue;
}

void LayerVolumeBase::SetBlankValue( float fBlank )
{
	m_fBlankValue = fBlank;
}

int LayerVolumeBase::GetBrushRadius()
{
	return m_nBrushRadius;
}

void LayerVolumeBase::SetBrushRadius( int nRadius )
{
	m_nBrushRadius = nRadius;
}
