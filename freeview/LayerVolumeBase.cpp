/**
 * @brief Layer data object for MRI volume.
 *
 */
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
 *
 */


#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "MyUtils.h"
#include "MyVTKUtils.h"
#include "BrushProperty.h"
#include "MainWindow.h"
#include "LivewireTool.h"
#include <stdlib.h>
#include <QFile>
#include <QDir>
#include <QDebug>
#include "vtkSimpleLabelEdgeFilter.h"
#include "vtkImageReslice.h"

LayerVolumeBase::LayerVolumeBase( QObject* parent ) : LayerEditable( parent )
{
  m_strTypeNames.push_back( "VolumeBase" );

  m_nActiveFrame = 0;
  m_propertyBrush = MainWindow::GetMainWindow()->GetBrushProperty();
  m_fFillValue = m_propertyBrush->GetFillValue();
  m_fBlankValue = m_propertyBrush->GetEraseValue();
  m_nBrushRadius = m_propertyBrush->GetBrushSize();
  m_b3DBrush = m_propertyBrush->Get3DBrush();
  m_livewire = new LivewireTool();
  m_imageData = NULL;
  m_imageDataRef = NULL;
  m_shiftBackgroundData = NULL;
  m_shiftForegroundData = NULL;
  connect(m_propertyBrush, SIGNAL(FillValueChanged(double)), this, SLOT(SetFillValue(double)));
  if (GetEndType() != "ROI")
    connect(m_propertyBrush, SIGNAL(EraseValueChanged(double)), this, SLOT(SetBlankValue(double)));
  connect(m_propertyBrush, SIGNAL(BrushSizeChanged(int)), this, SLOT(SetBrushRadius(int)));
  connect(m_propertyBrush, SIGNAL(Brush3DChanged(bool)), this, SLOT(Set3DBrush(bool)));
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

QVector<int> LayerVolumeBase::SetVoxelByIndex( int* n_in, int nPlane, bool bAdd, bool ignore_brush_size )
{
  QVector<int> indices;
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

  int nBrushSize = (ignore_brush_size? 1 : m_propertyBrush->GetBrushSize());
  int n[3], nsize[3] = { nBrushSize/2+1, nBrushSize/2+1, nBrushSize/2+1 };
  nsize[nPlane] = nBrushSize/2+1;
  int nActiveComp = GetActiveFrame();
  double* draw_range = bAdd ? m_propertyBrush->GetDrawRange(): m_propertyBrush->GetEraseRange();
  double* exclude_range = bAdd ? m_propertyBrush->GetExcludeRange() : m_propertyBrush->GetEraseExcludeRange();
  LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
  vtkImageData* ref = m_imageData;
  int nActiveCompRef = 0;
  if(!m_propertyBrush->Get3DBrush())
  {
      nsize[nPlane] = 1;
  }
  if ( ref_layer != NULL )
  {
    ref = ref_layer->GetImageData();
    nActiveCompRef = ref_layer->GetActiveFrame();
  }
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  char* ref_ptr = (char*)ref->GetScalarPointer();
  int* ref_dim = ref->GetDimensions();
  int ref_scalar_type = ref->GetScalarType();
  int ref_n_frames = ref->GetNumberOfScalarComponents();
  bool bNotROI = (GetEndType() != "ROI");
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[0], n[1], n[2], nActiveCompRef, ref_scalar_type );
          if (bAdd)
          {
            if ( ( bNotROI && m_propertyBrush->GetDrawRangeEnabled() &&
                   ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
                 ( bNotROI && m_propertyBrush->GetExcludeRangeEnabled() &&
                   ( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) ||
                 ( bNotROI && m_propertyBrush->GetDrawConnectedOnly() &&
                   ( !GetConnectedToOld( m_imageData, nActiveComp, n, nPlane ) ) ) )
            {
              ;
            }
            else
            {
              //   m_imageData->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fFillValue );
              MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, n[0], n[1], n[2], nActiveComp, scalar_type, m_fFillValue);
              indices << n[0] << n[1] << n[2];
              UpdateVoxelValueRange( m_fFillValue );
            }
          }
          else
          {
            if ( ( bNotROI && m_propertyBrush->GetEraseRangeEnabled() &&
                   ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
                 ( bNotROI && m_propertyBrush->GetEraseExcludeRangeEnabled() &&
                   ( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) )
            {
              ;
            }
            else
            {
              //  m_imageData->SetScalarComponentFromFloat( n[0], n[1], n[2], nActiveComp, m_fBlankValue );
              MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, n[0], n[1], n[2], nActiveComp, scalar_type, m_fBlankValue);
              indices << n[0] << n[1] << n[2];
            }
          }
        }
      }
    }
  }
  return indices;
}

bool LayerVolumeBase::CloneVoxelByIndex( int* n_in, int nPlane )
{
  int* nDim = m_imageData->GetDimensions();

  int nBrushSize = m_propertyBrush->GetBrushSize();
  int n[3], nsize[3] = { nBrushSize/2+1, nBrushSize/2+1, nBrushSize/2+1 };
  nsize[nPlane] = 1;
  int nActiveComp = GetActiveFrame();
  double* draw_range = m_propertyBrush->GetDrawRange();
  double* exclude_range = m_propertyBrush->GetExcludeRange();
  LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
  vtkImageData* ref = m_imageData;
  int nActiveCompRef = 0;
  if(!m_propertyBrush->Get3DBrush())
  {
      nsize[nPlane] = 1;
  }
  if ( ref_layer != NULL )
  {
    ref = ref_layer->GetImageData();
    nActiveCompRef = ref_layer->GetActiveFrame();
  }
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  char* ref_ptr = (char*)ref->GetScalarPointer();
  int* ref_dim = ref->GetDimensions();
  int ref_scalar_type = ref->GetScalarType();
  int ref_n_frames = ref->GetNumberOfScalarComponents();
  bool bNotROI = (GetEndType() != "ROI");
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[0], n[1], n[2], nActiveCompRef, ref_scalar_type );
          if ( ( bNotROI && m_propertyBrush->GetDrawRangeEnabled() &&
                 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
               ( bNotROI && m_propertyBrush->GetExcludeRangeEnabled() &&
                 ( fvalue >= exclude_range[0] && fvalue <= exclude_range[1] ) ) ||
               ( bNotROI && m_propertyBrush->GetDrawConnectedOnly() &&
                 ( !GetConnectedToOld( m_imageData, nActiveComp, n, nPlane ) ) ) )
          {
            ;
          }
          else
          {
            MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, n[0], n[1], n[2], nActiveComp, scalar_type, fvalue );
            UpdateVoxelValueRange( fvalue );
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
  char* ptr = (char*)img->GetScalarPointer();
  int scalar_type = img->GetScalarType();
  int n_frames = img->GetNumberOfScalarComponents();
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ptr, nDim, n_frames, n[0], n[1], n[2], nFrame, scalar_type );
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

void LayerVolumeBase::SetVoxelByRAS( double* ras, int nPlane, bool bAdd, bool ignore_brush_size )
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  QVector<int> list = SetVoxelByIndex( n, nPlane, bAdd, ignore_brush_size );
  if ( !list.isEmpty() )
  {
    SetModified();
    emit ActorUpdated();
    emit BaseVoxelEdited(list, bAdd);
  }
  else
  {
    // PopUndo();  // pop the previously saved undo step
  }
}

void LayerVolumeBase::SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd, bool ignore_brush_size )
{
  int n1[3], n2[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = ( int )( ( ras1[i] - origin[i] ) / voxel_size[i] + 0.5 );
    n2[i] = ( int )( ( ras2[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  QVector<int> list = SetVoxelByIndex( n1, n2, nPlane, bAdd, ignore_brush_size );
  if ( !list.isEmpty() )
  {
    SetModified();
    emit ActorUpdated();
    emit BaseVoxelEdited(list, bAdd);
  }
  else
  {
    // PopUndo();
  }
}

QVector<int> LayerVolumeBase::SetVoxelByIndex( int* n1, int* n2, int nPlane, bool bAdd, bool ignore_brush_size )
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
  QVector<int> list;
  list = SetVoxelByIndex( n1, nPlane, bAdd, ignore_brush_size );
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
      QVector<int> list1 = SetVoxelByIndex( n, nPlane, bAdd, ignore_brush_size );
      if (!list1.isEmpty())
        list << list1;
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
      QVector<int> list1 = SetVoxelByIndex( n, nPlane, bAdd, ignore_brush_size );
      if (!list1.isEmpty())
        list << list1;
    }
  }
  return list;
}

void LayerVolumeBase::CloneVoxelByRAS( double* ras, int nPlane )
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  if ( CloneVoxelByIndex( n, nPlane ) )
  {
    SetModified();
    emit ActorUpdated();
  }
  else
  {
    // PopUndo();  // pop the previously saved undo step
  }
}

void LayerVolumeBase::CloneVoxelByRAS( double* ras1, double* ras2, int nPlane )
{
  int n1[3], n2[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = ( int )( ( ras1[i] - origin[i] ) / voxel_size[i] + 0.5 );
    n2[i] = ( int )( ( ras2[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  if ( CloneVoxelByIndex( n1, n2, nPlane ) )
  {
    SetModified();
    emit ActorUpdated();
  }
  else
  {
    // PopUndo();
  }
}

bool LayerVolumeBase::CloneVoxelByIndex( int* n1, int* n2, int nPlane)
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
  bool bChanged = CloneVoxelByIndex( n1, nPlane );
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
      bChanged = CloneVoxelByIndex( n, nPlane );
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
      bChanged = CloneVoxelByIndex( n, nPlane );
    }
  }
  return true;
}

bool LayerVolumeBase::BorderFillByRAS(double *ras, int nPlane, bool b3D)
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  QVector<int> list = BorderFillByRAS(n, nPlane);
  if (!list.isEmpty())
  {
    SetModified();
    emit ActorUpdated();
    emit BaseVoxelEdited(list, true);
    return true;
  }
  else
    return false;
}

QVector<int> LayerVolumeBase::BorderFillByRAS(int *n, int nPlane)
{
  QVector<int> voxel_list;
  int* nDim = m_imageData->GetDimensions();
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int nx = 0, ny = 0;
  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  reslice->SetOutputDimensionality(2);
  double slicePos[3];
  for (int i = 0; i < 3; i++)
  {
    slicePos[i] = n[i]*voxel_size[i]+origin[i];
  }
  switch ( nPlane )
  {
  case 0:
    nx = nDim[1];
    ny = nDim[2];
    reslice->SetResliceAxesDirectionCosines( 0, 1, 0,
                                                 0, 0, 1,
                                                 1, 0, 0 );
    reslice->SetResliceAxesOrigin( slicePos[0], 0, 0  );
    break;
  case 1:
    nx = nDim[0];
    ny = nDim[2];
    reslice->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 0, 1,
                                                 0, 1, 0 );
    reslice->SetResliceAxesOrigin( 0, slicePos[1], 0 );
    break;
  case 2:
    nx = nDim[0];
    ny = nDim[1];
    reslice->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );
    reslice->SetResliceAxesOrigin( 0, 0, slicePos[2] );
    break;
  }

  LayerVolumeBase* ref_layer = m_propertyBrush->GetReferenceLayer();
  vtkImageData* ref = m_imageData;
  int nActiveCompRef = 0;
  if ( ref_layer != NULL )
  {
    ref = ref_layer->GetImageData();
    nActiveCompRef = ref_layer->GetActiveFrame();
  }
  int nActiveComp = this->GetActiveFrame();
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  char* ref_ptr = (char*)ref->GetScalarPointer();
  int* ref_dim = ref->GetDimensions();
  int ref_scalar_type = ref->GetScalarType();
  int ref_n_frames = ref->GetNumberOfScalarComponents();
  float fVoxelValue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[0], n[1], n[2], 0, ref_scalar_type);
  vtkSmartPointer<vtkSimpleLabelEdgeFilter> filter = vtkSmartPointer<vtkSimpleLabelEdgeFilter>::New();
#if VTK_MAJOR_VERSION > 5
  reslice->SetInputData(ref);
#else
  reslice->SetInput(ref);
#endif
  filter->SetInputConnection(reslice->GetOutputPort());
  filter->Update();
  vtkSmartPointer<vtkImageData> outline_image = filter->GetOutput();
  ref_ptr = (char*)outline_image->GetScalarPointer();
  ref_dim = outline_image->GetDimensions();
  float fTolerance = 0;
  switch ( nPlane )
  {
  case 0:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if (fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, j, 0, nActiveCompRef, scalar_type ) - fVoxelValue) <= fTolerance)
        {
          MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, n[nPlane], i, j, nActiveComp, scalar_type, m_fFillValue);
          voxel_list << n[nPlane] << i << j;
        }
      }
    }
    break;
  case 1:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if (fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, j, 0, nActiveCompRef, ref_scalar_type ) - fVoxelValue ) <= fTolerance)
        {
          MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, i, n[nPlane], j, nActiveComp, scalar_type, m_fFillValue);
          voxel_list << i << n[nPlane] << j;
        }
      }
    }
    break;
  case 2:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if (fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, j, 0, nActiveCompRef, ref_scalar_type ) - fVoxelValue ) <= fTolerance)
        {
          MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, i, j, n[nPlane], nActiveComp, scalar_type, m_fFillValue );
          voxel_list << i << j << n[nPlane];
        }
      }
    }
    break;
  }
  return voxel_list;
}

bool LayerVolumeBase::FloodFillByRAS( double* ras, int nPlane, bool bAdd, bool b3D, char* mask_out, bool ignore_exclusion )
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
    QVector<int> list = FloodFillByIndex( n, nPlane, bAdd, true, mask_out, ignore_exclusion );
    if ( !list.isEmpty() )
    {
      if ( !mask_out )
      {
        SetModified();
      }
      emit ActorUpdated();
      emit BaseVoxelEdited(list, bAdd);
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    QVector<int> list_all;
    int n0[3] = { n[0], n[1], n[2]};
    for (int i = n0[nPlane]; i < dim[nPlane]; i++)
    {
      n[nPlane] = i;
      QVector<int> list = FloodFillByIndex( n, nPlane, bAdd, false);
      if (list.isEmpty())
        break;
      else
        list_all << list;
    }
    for (int i = n0[nPlane]-1; i >= 0; i--)
    {
      n[nPlane] = i;
      QVector<int> list = FloodFillByIndex( n, nPlane, bAdd, false);
      if (list.isEmpty())
        break;
      else
        list_all << list;
    }
    SetModified();
    emit ActorUpdated();
    emit BaseVoxelEdited(list_all, bAdd);
    return true;
  }
}

// when mask_out is not null, do not fill the actual image data. instead, fill the mask_out buffer
QVector<int> LayerVolumeBase::FloodFillByIndex( int* n, int nPlane, bool bAdd, bool ignore_overflow, char* mask_out, bool ignore_exclusion )
{
  QVector<int> voxel_list;
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
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  char* ref_ptr = (char*)ref->GetScalarPointer();
  int* ref_dim = ref->GetDimensions();
  int ref_scalar_type = ref->GetScalarType();
  int ref_n_frames = ref->GetNumberOfScalarComponents();
  float fVoxelValue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[0], n[1], n[2], 0, ref_scalar_type);
  double fRange[2];
  ref->GetScalarRange( fRange );
  double fTolerance = (fRange[1]-fRange[0]) * m_propertyBrush->GetBrushTolerance() / 100.0;   // tolerance is percentage
  bool bNotROI = (GetEndType() != "ROI");
  switch ( nPlane )
  {
  case 0:
    if ( ref == m_imageData )
    {
      for ( i = 0; i < nx; i++ )
      {
        for ( j = 0; j < ny; j++ )
        {
          mask[j][i] = ( fabs( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, n[nPlane], i, j, nActiveCompRef, scalar_type ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, n[nPlane], i, j, nActiveComp, scalar_type ) <= 0  &&
                         fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[nPlane], i, j, nActiveCompRef, scalar_type ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( fabs( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, n[nPlane], j, nActiveCompRef, scalar_type ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, n[nPlane], j, nActiveComp, scalar_type ) <= 0  &&
                         fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, n[nPlane], j, nActiveCompRef, ref_scalar_type ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( fabs( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, j, n[nPlane], nActiveCompRef, scalar_type ) - fVoxelValue ) <= fTolerance
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
          mask[j][i] = ( MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, j, n[nPlane], nActiveComp, scalar_type ) <= 0  &&
                         fabs( MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, j, n[nPlane], nActiveCompRef, ref_scalar_type ) - fVoxelValue ) <= fTolerance
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
      return voxel_list;
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, n[nPlane], i, j, nActiveCompRef, ref_scalar_type );
          if ( ( bNotROI && m_propertyBrush->GetDrawRangeEnabled() &&
                 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
               ( bNotROI && !ignore_exclusion && m_propertyBrush->GetExcludeRangeEnabled() &&
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
              MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, n[nPlane], i, j, nActiveComp, scalar_type, bAdd ? m_fFillValue : m_fBlankValue );
            }
            voxel_list << n[nPlane] << i << j;
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, n[nPlane], j, nActiveCompRef, ref_scalar_type );
          if ( ( bNotROI && m_propertyBrush->GetDrawRangeEnabled() &&
                 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
               ( bNotROI && m_propertyBrush->GetExcludeRangeEnabled() &&
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
              MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, i, n[nPlane], j, nActiveComp, scalar_type, bAdd ? m_fFillValue : m_fBlankValue );
            }
            voxel_list << i << n[nPlane] << j;
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
          double fvalue = MyVTKUtils::GetImageDataComponent(ref_ptr, ref_dim, ref_n_frames, i, j, n[nPlane], nActiveCompRef, ref_scalar_type );
          if ( ( bNotROI && m_propertyBrush->GetDrawRangeEnabled() &&
                 ( fvalue < draw_range[0] || fvalue > draw_range[1] ) ) ||
               ( bNotROI && m_propertyBrush->GetExcludeRangeEnabled() &&
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
              MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, i, j, n[nPlane], nActiveComp, scalar_type, bAdd ? m_fFillValue : m_fBlankValue );
            }
            voxel_list << i << j << n[nPlane];
          }
        }
        ncnt++;
      }
    }
    break;
  }

  MyUtils::FreeMatrix( mask, ny );

  return voxel_list;
}

void LayerVolumeBase::SetLiveWireByRAS( double* pt1, double* pt2, int nPlane )
{
  int n1[3];
  double* orig = m_imageData->GetOrigin();
  double* vxlsize = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = ( int )( ( pt1[i] - orig[i] ) / vxlsize[i] );
    //    n2[i] = ( int )( ( pt2[i] - orig[i] ) / vxlsize[i] );
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
  QVector<int> list;
  for ( int i = 0; i < pts->GetNumberOfPoints(); i++ )
  {
    double* p = pts->GetPoint( i );
    n[0] = (int)( ( p[0] - orig[0] ) / vxlsize[0] + 0.5 );
    n[1] = (int)( ( p[1] - orig[1] ) / vxlsize[1] + 0.5);
    n[2] = (int)( ( p[2] - orig[2] ) / vxlsize[2] + 0.5 );

    list << SetVoxelByIndex( n, nPlane, true );
  }

  SetModified();
  emit ActorUpdated();
  emit BaseVoxelEdited(list, true);
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

  int n1[3];
  double* orig = image->GetOrigin();
  double* vxlsize = image->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = ( int )( ( pt1[i] - orig[i] ) / vxlsize[i] );
    //    n2[i] = ( int )( ( pt2[i] - orig[i] ) / vxlsize[i] );
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
  Q_UNUSED(n1);
  Q_UNUSED(n2);
  Q_UNUSED(nPlane);
  //  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
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
    SaveBufferItem( item2, item.plane, item.slice, item.frame );
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
    SaveBufferItem( item2, item.plane, item.slice, item.frame );
    m_bufferUndo.push_back( item2 );

    LoadBufferItem( item );
    item.Clear();

    SetModified();
    emit ActorUpdated();
    emit Modified();
  }
}

void LayerVolumeBase::SaveForUndo( int nPlane, bool bAllFrames )
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
  SaveBufferItem( item, nPlane, nSlice, bAllFrames? -1 : GetActiveFrame() );
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
  SaveBufferItem( m_bufferClipboard, nPlane, nSlice, GetActiveFrame() );
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
  long long nsize = ((long long)dim[0])*dim[1]*dim[2];
  char* mask = new char[nsize];
  memset( mask, 0, nsize );

  if ( FloodFillByRAS( ras, nPlane, true, false, mask, true ) )
  {
    m_bufferClipboard.Clear();

    SaveBufferItem( m_bufferClipboard, nPlane, nSlice[nPlane], GetActiveFrame(), mask );
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

void LayerVolumeBase::SaveBufferItem( UndoRedoBufferItem& item, int nPlane, int nSlice, int nFrame, const char* mask )
{
  item.plane = nPlane;
  item.slice = nSlice;
  item.frame = nFrame;
  if ( nPlane >= 0 )
  {
    int nDim[3], nStart[3] = { 0, 0, 0 };
    m_imageData->GetDimensions( nDim );
    nDim[nPlane] = 1;
    nStart[nPlane] = nSlice;
    size_t nSize = ((size_t)nDim[0])*nDim[1]*nDim[2]*m_imageData->GetScalarSize();
    item.data = new char[nSize];
    memset( item.data, 0, nSize );
    long long n = 0;
    char* ptr = (char*)m_imageData->GetScalarPointer();
    int scalar_size = m_imageData->GetScalarSize();
    int n_frames = m_imageData->GetNumberOfScalarComponents();
    int nOrigDim[3];
    m_imageData->GetDimensions( nOrigDim );
    for ( size_t i = nStart[0]; i < (size_t)nStart[0] + nDim[0]; i++ )
    {
      for ( size_t j = nStart[1]; j < (size_t)nStart[1] + nDim[1]; j++ )
      {
        for ( size_t k = nStart[2]; k < (size_t)nStart[2] + nDim[2]; k++ )
        {
          if ( !mask || mask[n] > 0 )
          {
            memcpy( item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * scalar_size,
                ptr + ((k*nOrigDim[0]*nOrigDim[1] + j*nOrigDim[0] + i) * n_frames + nFrame) * scalar_size,
                scalar_size );
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
    long long nSize = ((long long)nDim[0])*nDim[1]*nDim[2]*m_imageData->GetScalarSize()*(nFrame >= 0 ? 1 : m_imageData->GetNumberOfScalarComponents());
    QFile file(this->GenerateCacheFileName());
    if (!file.open(QIODevice::WriteOnly) || file.write((char*)m_imageData->GetScalarPointer() + (nFrame >= 0 ? nSize*nFrame : 0), nSize) != nSize)
    {
      cerr << "Could not write undo cache to disk" << endl;
      return;
    }
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
    char* ptr = (char*)m_imageData->GetScalarPointer();
    int scalar_size = m_imageData->GetScalarSize();
    int scalar_type = m_imageData->GetScalarType();
    int n_frames = m_imageData->GetNumberOfScalarComponents();
    int nOrigDim[3];
    m_imageData->GetDimensions( nOrigDim );
    for ( size_t i = nStart[0]; i < (size_t)nStart[0] + nDim[0]; i++ )
    {
      for ( size_t j = nStart[1]; j < (size_t)nStart[1] + nDim[1]; j++ )
      {
        for ( size_t k = nStart[2]; k < (size_t)nStart[2] + nDim[2]; k++ )
        {
          double dValue = MyVTKUtils::GetImageDataComponent(ptr, nOrigDim, n_frames, i, j, k, item.frame, scalar_type);
          memcpy( ptr + ((k*nOrigDim[0]*nOrigDim[1] + j*nOrigDim[0] + i)*n_frames + item.frame) * scalar_size,
              item.data + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * scalar_size,
              scalar_size );
          if ( bIgnoreZeros )
          {
            if ( MyVTKUtils::GetImageDataComponent(ptr, nOrigDim, n_frames, i, j, k, item.frame, scalar_type) == 0 )
            {
              MyVTKUtils::SetImageDataComponent(ptr, nOrigDim, n_frames, i, j, k, item.frame, scalar_type, dValue);
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
    int nSize = nDim[0]*nDim[1]*nDim[2]*m_imageData->GetScalarSize()*(item.frame >= 0 ? 1 : m_imageData->GetNumberOfScalarComponents());;
    char* p = (char*)m_imageData->GetScalarPointer() + (item.frame >= 0 ? nSize*item.frame : 0);
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

double LayerVolumeBase::GetFillValue()
{
  return m_fFillValue;
}

void LayerVolumeBase::SetFillValue( double fFill )
{
  if (m_fFillValue != fFill)
  {
    m_fFillValue = fFill;
    emit FillValueChanged( fFill );
  }
}

double LayerVolumeBase::GetBlankValue()
{
  return m_fBlankValue;
}

void LayerVolumeBase::SetBlankValue( double fBlank )
{
  if (m_fBlankValue != fBlank)
  {
    m_fBlankValue = fBlank;
    emit EraseValueChanged(fBlank);
  }
}

int LayerVolumeBase::GetBrushRadius()
{
  return m_nBrushRadius;
}

void LayerVolumeBase::SetBrushRadius( int nRadius )
{
  if (m_nBrushRadius != nRadius)
  {
    m_nBrushRadius = nRadius;
    emit BrushRadiusChanged(nRadius);
  }
}

void LayerVolumeBase::Set3DBrush( bool b3DBrush )
{
    if (m_b3DBrush != b3DBrush)
    {
        m_b3DBrush = b3DBrush;
        emit Brush3DChanged(m_b3DBrush);
    }
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
  QString strg = QDir::tempPath() + "/freeview-cache-" + QString::number(rand());
  while (QFile::exists(strg))
  {
    strg = QDir::tempPath() + "/freeview-cache-" + QString::number(rand());
  }
#ifdef Q_CYGWIN_WIN
  strg = MyUtils::NormalizeCygwinPath(strg);
#endif
  return strg;
}

void LayerVolumeBase::ClearVoxels()
{
  int* dim = m_imageData->GetDimensions();
  memset(m_imageData->GetScalarPointer(), 0, ((size_t)dim[0])*dim[1]*dim[2]*m_imageData->GetScalarSize());
  m_imageData->Modified();
}

void LayerVolumeBase::ShiftVoxelsByRAS(double* ras, int nPlane)
{
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }

  ShiftVoxels(n, nPlane);
}

void LayerVolumeBase::ShiftVoxels(int *nOffset, int nPlane)
{
  int nDim[3], nStart[3] = { 0, 0, 0 };
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  m_imageData->GetDimensions( nDim );
  nDim[nPlane] = 1;
  nStart[nPlane] = ( int )( ( m_dSlicePosition[nPlane] - origin[nPlane] ) / voxel_size[nPlane] + 0.5 );
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_size = m_imageData->GetScalarSize();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  int nFrame = m_nActiveFrame;
  int nOrigDim[3];
  m_imageData->GetDimensions( nOrigDim );
  for ( size_t i = nStart[0]; i < (size_t)nStart[0] + nDim[0]; i++ )
  {
    for ( size_t j = nStart[1]; j < (size_t)nStart[1] + nDim[1]; j++ )
    {
      for ( size_t k = nStart[2]; k < (size_t)nStart[2] + nDim[2]; k++ )
      {
        int ii = i-nStart[0]-nOffset[0];
        int jj = j-nStart[1]-nOffset[1];
        int kk = k-nStart[2]-nOffset[2];
        double val = 0;
        if (ii >= 0 && ii < nDim[0] && jj >= 0 && jj < nDim[1] && kk >= 0 && kk < nDim[2])
          val = MyVTKUtils::GetImageDataComponent(m_shiftForegroundData, nDim, 1, ii, jj, kk, 0, scalar_type);
        if (val > 0)
          MyVTKUtils::SetImageDataComponent(ptr, nOrigDim, n_frames, i, j, k, nFrame, scalar_type, m_fFillValue);
        else
          memcpy( ptr + ((k*nOrigDim[0]*nOrigDim[1] + j*nOrigDim[0] + i) * n_frames + nFrame) * scalar_size,
            m_shiftBackgroundData + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * scalar_size,
            scalar_size );
      }
    }
  }
  SetModified();
  emit ActorUpdated();
}

void LayerVolumeBase::PrepareShifting(int nPlane)
{
  delete[] m_shiftBackgroundData;
  delete[] m_shiftForegroundData;
  int nDim[3], nStart[3] = { 0, 0, 0 };
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  m_imageData->GetDimensions( nDim );
  nDim[nPlane] = 1;
  nStart[nPlane] = ( int )( ( m_dSlicePosition[nPlane] - origin[nPlane] ) / voxel_size[nPlane] + 0.5 );
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_size = m_imageData->GetScalarSize();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  int nFrame = m_nActiveFrame;
  size_t nsize = ((size_t)nDim[0])*nDim[1]*nDim[2]*scalar_size;
  m_shiftBackgroundData = new char[nsize];
  m_shiftForegroundData = new char[nsize];
  memset(m_shiftBackgroundData, 0, nsize);
  memset(m_shiftForegroundData, 0, nsize);
  int nOrigDim[3];
  m_imageData->GetDimensions( nOrigDim );
  for ( size_t i = nStart[0]; i < (size_t)nStart[0] + nDim[0]; i++ )
  {
    for ( size_t j = nStart[1]; j < (size_t)nStart[1] + nDim[1]; j++ )
    {
      for ( size_t k = nStart[2]; k < (size_t)nStart[2] + nDim[2]; k++ )
      {
        double val = MyVTKUtils::GetImageDataComponent(ptr, nOrigDim, n_frames, i, j, k, nFrame, scalar_type);
        if ( val == m_fFillValue )
        {
          MyVTKUtils::SetImageDataComponent(m_shiftForegroundData, nDim, 1, i-nStart[0], j-nStart[1], k-nStart[2], 0, scalar_type, 1);
        }
        else
        {
          memcpy( m_shiftBackgroundData + ( (k-nStart[2])*nDim[1]*nDim[0] + (j-nStart[1])*nDim[0] + (i-nStart[0]) ) * scalar_size,
              ptr + ((k*nOrigDim[0]*nOrigDim[1] + j*nOrigDim[0] + i) * n_frames + nFrame) * scalar_size,
              scalar_size );
        }
      }
    }
  }
}

void LayerVolumeBase::DoneShifting()
{

}
