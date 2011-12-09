/**
 * @file  LayerVolumeBase.cpp
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/09 22:09:05 $
 *    $Revision: 1.24 $
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


#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "MyUtils.h"
#include "BrushProperty.h"
#include "MainWindow.h"
#include "LivewireTool.h"
#include <stdlib.h>
#include <QFile>
#include <QDir>
#include <QDebug>

LayerVolumeBase::LayerVolumeBase( QObject* parent ) : LayerEditable( parent )
{
  m_strTypeNames.push_back( "VolumeBase" );

  m_fFillValue = 1;
  m_fBlankValue = 0;
  m_nBrushRadius = 1;
  m_nActiveFrame = 0;
  m_propertyBrush = MainWindow::GetMainWindow()->GetBrushProperty();
  m_livewire = new LivewireTool();
  m_imageData = NULL;
  m_imageDataRef = NULL;
}

LayerVolumeBase::~LayerVolumeBase()
{
  for ( size_t i = 0; i < m_bufferUndo.size(); i++ )
  {
    m_bufferUndo[i].Clear();
  }

  m_bufferUndo.clear();

  for ( size_t i = 0; i < m_bufferRedo.size(); i++ )
  {
    m_bufferRedo[i].Clear();
  }

  m_bufferClipboard.Clear();

  m_bufferUndo.clear();

  delete m_livewire;
}

bool LayerVolumeBase::SetVoxelByIndex( int* n_in, int nPlane, bool bAdd )
{
  int* nDim = m_imageData->GetDimensions();
  /* for ( int i = 0; i < 3; i++ )
   {
    if ( n_in[i] < 0 || n_in[i] >= nDim[i] )
     return false;
   }
  */
// float* ptr = ( float* )m_imageData->GetScalarPointer( n );
// if ( !ptr )
//  return false;

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
            {
              ;
            }
            else
            {
              m_imageData->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fFillValue );
              UpdateVoxelValueRange( m_fFillValue );
            }
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
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}

void LayerVolumeBase::SetVoxelByRAS( double* ras, int nPlane, bool bAdd )
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  if ( SetVoxelByIndex( n, nPlane, bAdd ) )
  {
    SetModified();
    emit ActorUpdated();
  }
  else
  {
    // PopUndo();  // pop the previously saved undo step
  }
}

void LayerVolumeBase::SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd )
{
  int n1[3], n2[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = ( int )( ( ras1[i] - origin[i] ) / voxel_size[i] + 0.5 );
    n2[i] = ( int )( ( ras2[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  if ( SetVoxelByIndex( n1, n2, nPlane, bAdd ) )
  {
    SetModified();
    emit ActorUpdated();
  }
  else
  {
    // PopUndo();
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

bool LayerVolumeBase::FloodFillByRAS( double* ras, int nPlane, bool bAdd, bool b3D, char* mask_out )
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  if (!b3D)
  {
    if ( FloodFillByIndex( n, nPlane, bAdd, true, mask_out ) )
    {
      if ( !mask_out )
      {
        SetModified();
      }
      emit ActorUpdated();
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    int n0[3] = { n[0], n[1], n[2]};
    for (int i = n0[nPlane]; i < dim[nPlane]; i++)
    {
      n[nPlane] = i;
      if (!FloodFillByIndex( n, nPlane, bAdd, false))
        break;
    }
    for (int i = n0[nPlane]-1; i >= 0; i--)
    {
      n[nPlane] = i;
      if (!FloodFillByIndex( n, nPlane, bAdd, false))
        break;
    }
    SetModified();
    emit ActorUpdated();
    return true;
  }
}

// when mask_out is not null, do not fill the actual image data. instead, fill the mask_out buffer
bool LayerVolumeBase::FloodFillByIndex( int* n, int nPlane, bool bAdd, bool ignore_overflow, char* mask_out )
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
  double fRange[2];
  ref->GetScalarRange( fRange );
  double fTolerance = (fRange[1]-fRange[0]) * m_propertyBrush->GetBrushTolerance() / 100.0;   // tolerance is percentage
  switch ( nPlane )
  {
  case 0:
    if ( ref == m_imageData )
    {
      for ( i = 0; i < nx; i++ )
      {
        for ( j = 0; j < ny; j++ )
        {
          mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveCompRef ) - fVoxelValue ) <= fTolerance
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
                         fabs( ref->GetScalarComponentAsFloat( n[nPlane], i, j, nActiveCompRef ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveCompRef ) - fVoxelValue ) <= fTolerance
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
                         fabs( ref->GetScalarComponentAsFloat( i, n[nPlane], j, nActiveCompRef ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( fabs( ref->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveCompRef ) - fVoxelValue ) <= fTolerance
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
                         fabs( ref->GetScalarComponentAsFloat( i, j, n[nPlane], nActiveCompRef ) - fVoxelValue ) <= fTolerance
                         ? 1 : 0 );
        }
      }
    }
    break;
  }

  MyUtils::FloodFill( mask, x, y, 0, 0, nx-1, ny-1, 2, 0 );
  if (!ignore_overflow)
  {
    if (mask[0][0] == 2 && mask[ny-1][nx=1] == 2)
    {
      MyUtils::FreeMatrix( mask, ny );
      return false;
    }
  }
  int ncnt;
  switch ( nPlane )
  {
  case 0:
    ncnt = 0;
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
          {
            ;
          }
          else
          {
            if ( mask_out )
            {
              mask_out[ncnt] = 1;
            }
            else
            {
              m_imageData->SetScalarComponentFromFloat( n[nPlane], i, j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
            }
          }
        }
        ncnt++;
      }
    }
    break;
  case 1:
    ncnt = 0;
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
          {
            ;
          }
          else
          {
            if ( mask_out )
            {
              mask_out[ncnt] = 1;
            }
            else
            {
              m_imageData->SetScalarComponentFromFloat( i, n[nPlane], j, nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
            }
          }
        }
        ncnt++;
      }
    }
    break;
  case 2:
    ncnt = 0;
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
          {
            ;
          }
          else
          {
            if ( mask_out )
            {
              mask_out[ncnt] = 1;
            }
            else
            {
              m_imageData->SetScalarComponentFromFloat( i, j, n[nPlane], nActiveComp, bAdd ? m_fFillValue : m_fBlankValue );
            }
          }
        }
        ncnt++;
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
  {
    image = ref_layer->GetImageData();
  }
// else if ( m_imageDataRef.GetPointer() != NULL )
  else if ( m_imageDataRef != NULL )
  {
    image = m_imageDataRef;
  }

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
  emit ActorUpdated();
}

std::vector<double> LayerVolumeBase::GetLiveWirePointsByRAS( double* pt1, double* pt2, int nPlane )
{
  vtkImageData* image = m_imageData;
  LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
  if ( ref_layer != NULL )
  {
    image = ref_layer->GetImageData();
  }
// else if ( m_imageDataRef.GetPointer() != NULL )
  else if ( m_imageDataRef != NULL )
  {
    image = m_imageDataRef;
  }

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
// MyUtils::GetLivewirePoints( m_imageData, nPlane, n1[nPlane], n1, n2, pts );
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
    item.Clear();

    SetModified();
    emit ActorUpdated();
    emit Modified();
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
    item.Clear();

    SetModified();
    emit ActorUpdated();
    emit Modified();
  }
}

void LayerVolumeBase::SaveForUndo( int nPlane )
{
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int nSlice = 0;
  if ( nPlane >= 0 )
  {
    nSlice = ( int )( ( m_dSlicePosition[nPlane] - origin[nPlane] ) / voxel_size[nPlane] + 0.5 );
  }

  if ( (int)m_bufferUndo.size() >= m_nMaxUndoSteps )
  {
    m_bufferUndo[0].Clear();
    m_bufferUndo.erase( m_bufferUndo.begin() );
  }

  UndoRedoBufferItem item;
  SaveBufferItem( item, nPlane, nSlice );
  m_bufferUndo.push_back( item );

  // clear redo buffer
  for ( size_t i = 0; i < m_bufferRedo.size(); i++ )
  {
    m_bufferRedo[i].Clear();
  }
  m_bufferRedo.clear();
}

bool LayerVolumeBase::IsValidToPaste( int nPlane )
{
  return ( m_bufferClipboard.data != NULL && m_bufferClipboard.plane == nPlane );
}

void LayerVolumeBase::Copy( int nPlane )
{
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int nSlice = ( int )( ( m_dSlicePosition[nPlane] - origin[nPlane] ) / voxel_size[nPlane] + 0.5 );
  m_bufferClipboard.Clear();
  SaveBufferItem( m_bufferClipboard, nPlane, nSlice );
}

bool LayerVolumeBase::CopyStructure( int nPlane, double* ras )
{
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int dim[3];
  m_imageData->GetDimensions( dim );
  int nSlice[3];
  for ( int i = 0; i < 3; i++ )
  {
    nSlice[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
    if ( nSlice[i] < 0 || nSlice[i] >= dim[i] )
    {
      return false;
    }
  }
  if ( m_imageData->GetScalarComponentAsDouble( nSlice[0], nSlice[1], nSlice[2], 0 ) == 0 )
  {
    return false;
  }

  dim[nPlane] = 1;
  char* mask = new char[dim[0]*dim[1]*dim[2]];
  memset( mask, 0, dim[0]*dim[1]*dim[2] );

  if ( FloodFillByRAS( ras, nPlane, true, mask ) )
  {
    m_bufferClipboard.Clear();

    SaveBufferItem( m_bufferClipboard, nPlane, nSlice[nPlane], mask );
    return true;
  }
  else
  {
    return false;
  }
}

void LayerVolumeBase::Paste( int nPlane )
{
  SaveForUndo( nPlane );

  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int nSlice = ( int )( ( m_dSlicePosition[nPlane] - origin[nPlane] ) / voxel_size[nPlane] + 0.5 );
  m_bufferClipboard.slice = nSlice;
  LoadBufferItem( m_bufferClipboard, true );   // ignore zeros

  SetModified();
  emit ActorUpdated();
}

void LayerVolumeBase::SaveBufferItem( UndoRedoBufferItem& item, int nPlane, int nSlice, const char* mask )
{
  item.plane = nPlane;
  item.slice = nSlice;
  if ( nPlane >= 0 )
  {
    int nDim[3], nStart[3] = { 0, 0, 0 };
    m_imageData->GetDimensions( nDim );
    nDim[nPlane] = 1;
    nStart[nPlane] = nSlice;
    int nSize = nDim[0]*nDim[1]*nDim[2]*m_imageData->GetScalarSize();
    item.data = new char[nSize];
    memset( item.data, 0, nSize );
    int n = 0;
    for ( int i = nStart[0]; i < nStart[0] + nDim[0]; i++ )
    {
      for ( int j = nStart[1]; j < nStart[1] + nDim[1]; j++ )
      {
        for ( int k = nStart[2]; k < nStart[2] + nDim[2]; k++ )
        {
          if ( !mask || mask[n] > 0 )
          {
            memcpy( item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_imageData->GetScalarSize(),
                    m_imageData->GetScalarPointer( i, j, k ),
                    m_imageData->GetScalarSize() );
          }
          n++;
        }
      }
    }
  }
  else    // save whole volume in cache file
  {
    int nDim[3];
    m_imageData->GetDimensions( nDim );
    int nSize = nDim[0]*nDim[1]*nDim[2]*m_imageData->GetScalarSize();
    QFile file(this->GenerateCacheFileName());
    if (!file.open(QIODevice::WriteOnly))
    {
      return;
    }
    file.write((char*)m_imageData->GetScalarPointer() + nSize*m_nActiveFrame, nSize);
    file.close();
    item.cache_filename = file.fileName();
    if (this->IsTypeOf("MRI"))
    {
      LayerMRI* mri = qobject_cast<LayerMRI*>(this);
      item.mri_settings = mri->GetProperty()->GetSettings();
    }
  }
}

void LayerVolumeBase::LoadBufferItem( UndoRedoBufferItem& item, bool bIgnoreZeros )
{
  if (item.plane >= 0)
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
          double dValue = m_imageData->GetScalarComponentAsDouble( i, j, k, 0 );
          memcpy( m_imageData->GetScalarPointer( i, j, k ),
                  item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * m_imageData->GetScalarSize(),
                  m_imageData->GetScalarSize() );
          if ( bIgnoreZeros )
          {
            if ( m_imageData->GetScalarComponentAsDouble( i, j, k, 0 ) == 0 )
            {
              m_imageData->SetScalarComponentFromDouble( i, j, k, 0, dValue );
            }
          }
        }
      }
    }
  }
  else if (!item.cache_filename.isEmpty())
  {
    QFile file(item.cache_filename);
    if (!file.open(QIODevice::ReadOnly))
    {
      return;
    }
    int nDim[3];
    m_imageData->GetDimensions( nDim );
    int nSize = nDim[0]*nDim[1]*nDim[2]*m_imageData->GetScalarSize();
    char* p = (char*)m_imageData->GetScalarPointer() + nSize*m_nActiveFrame;
    file.read(p, nSize);
    file.close();
    if (this->IsTypeOf("MRI"))
    {
      LayerMRI* mri = qobject_cast<LayerMRI*>(this);
      mri->GetProperty()->RestoreSettings(item.mri_settings);
    }
    m_imageData->Modified();
  }
}

float LayerVolumeBase::GetFillValue()
{
  return m_fFillValue;
}

void LayerVolumeBase::SetFillValue( float fFill )
{
  m_fFillValue = fFill;
  emit FillValueChanged( fFill );
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

double LayerVolumeBase::GetMinimumVoxelSize()
{
  double* vs = m_imageData->GetSpacing();
  if ( vs[0] < vs[1] && vs[0] < vs[2] )
  {
    return vs[0];
  }
  else if ( vs[1] < vs[2] )
  {
    return vs[1];
  }
  else
  {
    return vs[2];
  }
}

void LayerVolumeBase::GetDisplayBounds( double* bounds )
{
  m_imageData->GetBounds( bounds );
}

QString LayerVolumeBase::GenerateCacheFileName()
{
  QString strg = QDir::tempPath() + "/freeview-cache-" + QString::number(qrand());
  while (QFile::exists(strg))
  {
    strg = QDir::tempPath() + "/freeview-cache-" + QString::number(qrand());
  }
#ifdef Q_CYGWIN_WIN
  strg = MyUtils::NormalizeCygwinPath(strg);
#endif
  return strg;
}
