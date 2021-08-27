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

#include "LayerROI.h"
#include "vtkRenderer.h"
#include "vtkImageReslice.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageFlip.h"
#include "vtkTransform.h"
#include "vtkPlaneSource.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTexture.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkImageChangeInformation.h"
#include "vtkImageMapper3D.h"
#include "vtkImageDilateErode3D.h"
#include "vtkImageOpenClose3D.h"
#include "LayerPropertyROI.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "FSLabel.h"
#include "FSSurface.h"
#include <stdlib.h>
#include <QDebug>
#include "LayerSurface.h"
#include "FSVolume.h"
#include "vtkImageThreshold.h"
#include "MyVTKUtils.h"

LayerROI::LayerROI( LayerMRI* layerMRI, QObject* parent ) : LayerVolumeBase( parent )
{
  m_strTypeNames << "Supplement" << "ROI";
  m_sPrimaryType = "ROI";
  m_nVertexCache = NULL;

  m_label = new FSLabel( this, layerMRI->GetSourceVolume() );
  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = 0;
    m_sliceActor2D[i] = vtkImageActor::New();
    m_sliceActor3D[i] = vtkImageActor::New();
    m_sliceActor2D[i]->InterpolateOff();
    m_sliceActor3D[i]->InterpolateOff();
#if VTK_MAJOR_VERSION > 5
    m_sliceActor2D[i]->ForceOpaqueOn();
    m_sliceActor3D[i]->ForceOpaqueOn();
#endif
  }

  mProperty = new LayerPropertyROI( this );

  m_layerMappedSurface = NULL;
  m_layerSource = layerMRI;
  m_imageDataRef = layerMRI->GetImageData();
  if ( m_layerSource )
  {
    SetWorldOrigin( m_layerSource->GetWorldOrigin() );
    SetWorldVoxelSize( m_layerSource->GetWorldVoxelSize() );
    SetWorldSize( m_layerSource->GetWorldSize() );

    m_imageData = vtkSmartPointer<vtkImageData>::New();
    // m_imageData->DeepCopy( m_layerSource->GetRASVolume() );

    m_imageData->SetOrigin( GetWorldOrigin() );
    m_imageData->SetSpacing( GetWorldVoxelSize() );
    m_imageData->SetDimensions( ( int )( m_dWorldSize[0] / m_dWorldVoxelSize[0] + 0.5 ),
        ( int )( m_dWorldSize[1] / m_dWorldVoxelSize[1] + 0.5 ),
        ( int )( m_dWorldSize[2] / m_dWorldVoxelSize[2] + 0.5 ) );
#if VTK_MAJOR_VERSION > 5
    m_imageData->AllocateScalars(VTK_FLOAT, 1);
#else
    m_imageData->SetScalarTypeToFloat();
    m_imageData->SetNumberOfScalarComponents(1);
    m_imageData->AllocateScalars();
#endif
    float* ptr = (float*)m_imageData->GetScalarPointer();
    int* dim = m_imageData->GetDimensions();
    size_t nsize = ((size_t)dim[0])*dim[1]*dim[2];
    for (size_t i = 0; i < nsize; i++)
    {
      ptr[i] = -1;
    }
    m_fBlankValue = -1;
    InitializeActors();
  }

  connect( mProperty, SIGNAL(ColorMapChanged()), this, SLOT(UpdateColorMap()) );
  connect( mProperty, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()) );
  connect( mProperty, SIGNAL(ThresholdChanged(double)), this, SLOT(UpdateThreshold()));

  connect(this, SIGNAL(BaseVoxelEdited(QVector<int>,bool)), SLOT(OnBaseVoxelEdited(QVector<int>,bool)));
  UpdateProperties();
}

LayerROI::~LayerROI()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->Delete();
    m_sliceActor3D[i]->Delete();
  }
  if (m_nVertexCache)
    delete[] m_nVertexCache;
}

bool LayerROI::LoadROIFromFile( const QString& filename )
{
  if ( !m_label->LabelRead( filename.toLatin1().data() ) )
  {
    return false;
  }
  UpdateProperties();
  m_label->Initialize(m_layerSource->GetSourceVolume(), 0, 0);
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume(), GetProperty()->GetThreshold());

  m_sFilename = filename;

  return true;
}

void LayerROI::UpdateProperties()
{
  double range[2];
  m_label->GetStatsRange(range);
  if (range[0] == range[1])
    range[1] = range[0] + 1;
  GetProperty()->SetHeatscaleValues(range[0], range[1]);
  GetProperty()->SetValueRange(range);
  m_fBlankValue = range[0]-1;
  GetProperty()->UpdateLUTTable();
}

void LayerROI::InitializeActors()
{
  if ( m_layerSource == NULL )
  {
    return;
  }

  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
#if VTK_MAJOR_VERSION > 5
    mReslice[i]->SetInputData( m_imageData );
#else
    mReslice[i]->SetInput( m_imageData );
#endif
    //  mReslice[i]->SetOutputSpacing( sizeX, sizeY, sizeZ );
    mReslice[i]->BorderOff();

    // This sets us to extract slices.
    mReslice[i]->SetOutputDimensionality( 2 );

    // This will change depending what orienation we're in.
    mReslice[i]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );

    // This will change to select a different slice.
    mReslice[i]->SetResliceAxesOrigin( 0, 0, 0 );
    //
    // Image to colors using color table.
    //
    mColorMap[i] = vtkSmartPointer<vtkImageMapToColors>::New();
    mColorMap[i]->SetLookupTable( GetProperty()->GetLookupTable() );
    mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    mColorMap[i]->SetOutputFormatToRGBA();
    mColorMap[i]->PassAlphaToOutputOn();

    //
    // Prop in scene with plane mesh and texture.
    //
    m_sliceActor2D[i]->GetMapper()->SetInputConnection( mColorMap[i]->GetOutputPort() );
    m_sliceActor3D[i]->GetMapper()->SetInputConnection( mColorMap[i]->GetOutputPort() );

    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }

  this->UpdateOpacity();
  this->UpdateColorMap();
}


void LayerROI::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetOpacity( GetProperty()->GetOpacity() );
    m_sliceActor3D[i]->SetOpacity( GetProperty()->GetOpacity() );
  }

  emit ActorUpdated();
}

void LayerROI::UpdateColorMap ()
{
  for ( int i = 0; i < 3; i++ )
  {
    mColorMap[i]->SetLookupTable( GetProperty()->GetLookupTable() );
    //    m_sliceActor2D[i]->GetProperty()->SetColor(1, 0, 0);
  }

  emit ActorUpdated();
}

void LayerROI::UpdateThreshold()
{
  if (m_label)
  {
    m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume(), GetProperty()->GetThreshold() );
    emit ActorUpdated();
  }
}

void LayerROI::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerROI::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( bSliceVisibility == NULL || bSliceVisibility[i] )
    {
      renderer->AddViewProp( m_sliceActor3D[i] );
    }
  }
}

void LayerROI::OnSlicePositionChanged( int nPlane )
{
  if ( !m_layerSource )
  {
    return;
  }

  vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  switch ( nPlane )
  {
  case 0:
    m_sliceActor2D[0]->PokeMatrix( matrix );
    m_sliceActor2D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor2D[0]->RotateX( 90 );
    m_sliceActor2D[0]->RotateY( -90 );
    m_sliceActor3D[0]->PokeMatrix( matrix );
    m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor3D[0]->RotateX( 90 );
    m_sliceActor3D[0]->RotateY( -90 );

    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice[0]->SetResliceAxesDirectionCosines( 0, -1, 0,
                                                 0, 0, 1,
                                                 1, 0, 0 );
    mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
    mReslice[0]->Modified();
    break;
  case 1:
    m_sliceActor2D[1]->PokeMatrix( matrix );
    m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor2D[1]->RotateX( 90 );
    // m_sliceActor2D[1]->RotateY( 180 );
    m_sliceActor3D[1]->PokeMatrix( matrix );
    m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor3D[1]->RotateX( 90 );
    // m_sliceActor3D[1]->RotateY( 180 );

    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice[1]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 0, 1,
                                                 0, 1, 0 );
    mReslice[1]->SetResliceAxesOrigin( 0, m_dSlicePosition[1], 0 );
    mReslice[1]->Modified();
    break;
  case 2:
    m_sliceActor2D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
    // m_sliceActor2D[2]->RotateY( 180 );
    m_sliceActor3D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
    // m_sliceActor3D[2]->RotateY( 180 );

    mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );
    mReslice[2]->SetResliceAxesOrigin( 0, 0, m_dSlicePosition[2]  );
    mReslice[2]->Modified();
    break;
  }
}

void LayerROI::SetVisible( bool bVisible )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
  }
  LayerVolumeBase::SetVisible(bVisible);
}

bool LayerROI::IsVisible()
{
  return m_sliceActor2D[0]->GetVisibility() > 0;
}

void LayerROI::SetModified()
{
  mReslice[0]->Modified();
  mReslice[1]->Modified();
  mReslice[2]->Modified();

  LayerVolumeBase::SetModified();
}

bool LayerROI::SaveROI( )
{
  if ( m_sFilename.size() == 0 || m_imageData.GetPointer() == NULL )
  {
    return false;
  }

  //    if (!m_layerMappedSurface)
  //        m_label->UpdateLabelFromImage( m_imageData, m_layerSource->GetSourceVolume() );

  bool bSaved = m_label->LabelWrite( m_sFilename.toLatin1().data() );
  if ( !bSaved )
  {
    m_bModified = true;
  }

  return bSaved;
}

bool LayerROI::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( m_sliceActor2D[i] == prop || m_sliceActor3D[i] == prop )
    {
      return true;
    }
  }
  return false;
}

bool LayerROI::DoRotate( std::vector<RotationElement>& rotations )
{
  Q_UNUSED(rotations);
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume() );
  return true;
}

void LayerROI::DoRestore()
{
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume() );
}

void LayerROI::UpdateLabelData()
{
  if ( IsModified() )
  {
    m_label->UpdateLabelFromImage( m_imageData, m_layerSource->GetSourceVolume() );
  }
}

bool LayerROI::GetCentroidPosition(double *pos)
{
  //    UpdateLabelData();
  if (m_label->GetCentroidRASPosition(pos, m_layerSource->GetSourceVolume()))
  {
    m_layerSource->RASToTarget(pos, pos);
    return true;
  }
  else
    return false;
}

void LayerROI::GetStats(int nPlane, int *count_out, float *area_out,
                        LayerMRI *underlying_mri, double *mean_out, double *sd_out)
{
  if ( !m_imageData || nPlane < 0 || nPlane > 2 )
  {
    return;
  }

  int* dim = m_imageData->GetDimensions();
  double* origin = m_imageData->GetOrigin();
  double vs[3];
  m_imageData->GetSpacing( vs );
  float* ptr = (float*)m_imageData->GetScalarPointer();

  int cnt = 0;
  //  QVector<int> indices;
  double val_range[2];
  GetProperty()->GetValueRange(val_range);
  QVector<float> coords;
  for ( int i = 0; i < dim[0]; i++ )
  {
    for ( int j = 0; j < dim[1]; j++ )
    {
      for ( int k = 0; k < dim[2]; k++ )
      {
        if ( ptr[k*dim[0]*dim[1]+j*dim[0]+i] >= val_range[0] )
        {
          cnt++;
          //          indices << i << j << k;
          coords << i*vs[0]+origin[0] << j*vs[1]+origin[1] << k*vs[2]+origin[2];
        }
      }
    }
  }
  //  vs[nPlane] = 1.0;

  *count_out = cnt;
  *area_out = cnt*vs[0]*vs[1]*vs[2];

  if (underlying_mri)
    underlying_mri->GetVoxelStatsByTargetRAS(coords, mean_out, sd_out);
}

void LayerROI::SetMappedSurface(LayerSurface *s)
{
  if (m_layerMappedSurface)
    m_layerMappedSurface->RemoveMappedLabel(this);
  m_layerMappedSurface = s;
  if (s)
  {
    s->AddMappedLabel(this);
    m_label->Initialize(m_layerSource->GetSourceVolume(), s->GetSourceSurface(), s->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    OnUpdateLabelRequested();
    SetModified();
  }
  else
  {
    m_label->Initialize(m_layerSource->GetSourceVolume(), NULL, 0);
  }
}

void LayerROI::OnSurfaceDestroyed(QObject *obj)
{
  if (m_layerMappedSurface == obj)
  {
    m_layerMappedSurface = NULL;
    m_label->Initialize(m_layerSource->GetSourceVolume(), NULL, 0);
  }
}

void LayerROI::OnUpdateLabelRequested()
{
  //    UpdateLabelData();
  if (m_layerMappedSurface)
  {
    //        int coords = CURRENT_VERTICES;
    //        if (m_layerMappedSurface->IsInflated())
    //            coords = WHITE_VERTICES;
    //        m_label->FillUnassignedVertices(m_layerMappedSurface->GetSourceSurface(), m_layerSource->GetSourceVolume(), coords);
    m_layerMappedSurface->UpdateOverlay();
  }
}

void LayerROI::MapLabelColorData( unsigned char* colordata, int nVertexCount )
{
  double* rgbColor = GetProperty()->GetColor();
  double dThreshold = GetProperty()->GetThreshold();
  LABEL* label = m_label->GetRawLabel();
  for ( int i = 0; i < label->n_points; i++ )
  {
    int vno = label->lv[i].vno;
    if (vno < nVertexCount && vno >= 0 && !label->lv[i].deleted)
    {
      double opacity = 1;
      if (label->lv[i].stat >= dThreshold)
        opacity = 1;
      else
        opacity = 0;
      double rgb[4] = { rgbColor[0], rgbColor[1], rgbColor[2], 1 };
      if (GetProperty()->GetColorCode() == LayerPropertyROI::Heatscale)
      {
        GetProperty()->GetLookupTable()->GetColor(label->lv[i].stat, rgb);
      }
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + rgb[0] * 255 * opacity );
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + rgb[1] * 255 * opacity );
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + rgb[2] * 255 * opacity );
    }
  }
}

void LayerROI::OnBaseVoxelEdited(const QVector<int>& voxel_list, bool bAdd)
{
  if (true)
  {
    int total_cnt = 0;
    for (int i = 0; i < voxel_list.size(); i+=3)
    {
      int n[3] = {voxel_list[i], voxel_list[i+1], voxel_list[i+2]};
      m_layerSource->TargetIndexToOriginalIndex(n, n);
      if (!m_nVertexCache && m_layerMappedSurface)
        m_nVertexCache = new int[m_layerMappedSurface->GetNumberOfVertices()];

      int cnt = 0;
      int coords = 0;
      if (m_layerMappedSurface)
        coords = m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES;
      m_label->EditVoxel(n[0], n[1], n[2], coords, bAdd, m_nVertexCache, m_nVertexCache?(&cnt):NULL);

      total_cnt += cnt;
    }
    if (total_cnt > 0 && m_layerMappedSurface)
    {
      //    qDebug() << "# of edited vertices:" << total_cnt;
      m_layerMappedSurface->UpdateOverlay(true, true);
    }
  }
}

void LayerROI::EditVertex(int nvo, bool bAdd)
{
  QVector<int> list;
  list << nvo;
  EditVertex(list, bAdd);
}

void LayerROI::OnLabelDataUpdated()
{
  if (m_layerMappedSurface)
  {
    if (m_label->UpdateStatsRange(0))
      UpdateProperties();
    m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume(), GetProperty()->GetThreshold());
    m_layerMappedSurface->UpdateOverlay(true, true);
  }
  else
  {
    //  m_label->UpdateLabelFromImage(m_imageData, m_layerSource->GetSourceVolume());
  }
  SetModified();
  emit ActorUpdated();
}

void LayerROI::EditVertex(const QVector<int> list_nvo_in, bool bAdd)
{
  if (m_layerMappedSurface)
  {
    QVector<int> list_nvo;
    if (m_nBrushRadius == 1)
      list_nvo = list_nvo_in;
    else
    {
      m_layerMappedSurface->SetNeighborhoodSize(qMin(3, m_nBrushRadius-1));
      for (int i = 0; i < list_nvo_in.size(); i++)
        list_nvo << list_nvo_in[i] << m_layerMappedSurface->GetVertexNeighbors(list_nvo_in[i]);
    }

    int ret = -1;
    int coords = m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES;
    for (int i = 0; i < list_nvo.size(); i++)
    {
      int nvo = list_nvo[i];
      if (bAdd)
        ret = qMax(ret, ::LabelAddVertex(m_label->GetRawLabel(), nvo, coords));
      else
        ret = qMax(ret, ::LabelDeleteVertex(m_label->GetRawLabel(), nvo, coords));
    }
    if (ret >= 0)
    {
      OnLabelDataUpdated();
    }
  }
}

void LayerROI::UpdateFilteredImage(vtkImageData* mask_before, vtkImageData* mask_after)
{
  int scalar_type = m_imageData->GetScalarType();
  char* ptr1 = (char*)mask_before->GetScalarPointer();
  char* ptr2 = (char*)mask_after->GetScalarPointer();
  char* ptr_img = (char*)m_imageData->GetScalarPointer();
  QVector<int> listAdd, listRemove;
  int* dim = m_imageData->GetDimensions();
  for (size_t i = 0; i < dim[0]; i++)
  {
    for (size_t j = 0; j < dim[1]; j++)
    {
      for (size_t k = 0; k < dim[2]; k++)
      {
        double val1 = MyVTKUtils::GetImageDataComponent(ptr1, dim, 1, i, j, k, 0, scalar_type);
        double val2 = MyVTKUtils::GetImageDataComponent(ptr2, dim, 1, i, j, k, 0, scalar_type);
        if (val1 == 0 && val2 >= 1)
        {
          listAdd << i << j << k;
          MyVTKUtils::SetImageDataComponent(ptr_img, dim, 1, i, j, k, 0, scalar_type, m_fFillValue);
        }
        else if (val1 == 2 && val2 <= 1)
        {
          listRemove << i << j << k;
          MyVTKUtils::SetImageDataComponent(ptr_img, dim, 1, i, j, k, 0, scalar_type, m_fBlankValue);
        }
      }
    }
  }
  if (!listAdd.isEmpty())
    OnBaseVoxelEdited(listAdd, true);
  if (!listRemove.isEmpty())
    OnBaseVoxelEdited(listRemove, true);
  m_imageData->Modified();
}

vtkSmartPointer<vtkImageData> LayerROI::GetThresholdedMaskImage()
{
  double range[2];
  m_label->GetStatsRange(range);
  vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
  threshold->ThresholdByUpper(range[0]);
  threshold->SetInValue(2);
  threshold->SetOutValue(0);
  threshold->ReplaceInOn();
  threshold->ReplaceOutOn();
  threshold->SetOutputScalarType(m_imageData->GetScalarType());
#if VTK_MAJOR_VERSION > 5
  threshold->SetInputData(m_imageData);
#else
  threshold->SetInput(m_imageData);
#endif
  threshold->Update();
  vtkSmartPointer<vtkImageData> out = threshold->GetOutput();
  return out;
}

void LayerROI::Dilate(int nTimes)
{
  SaveForUndo();
  if (m_layerMappedSurface)
  {
    ::LabelDilate(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes,
                  m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
  }
  else
  {
    vtkSmartPointer<vtkImageData> data1 = GetThresholdedMaskImage();
    vtkSmartPointer<vtkImageData> data2 = data1;
    for (int i = 0; i < nTimes; i++)
    {
      vtkSmartPointer<vtkImageDilateErode3D> filter = vtkSmartPointer<vtkImageDilateErode3D>::New();
      filter->SetKernelSize(2,2,2);
      filter->SetDilateValue(2);
      filter->SetErodeValue(0);
#if VTK_MAJOR_VERSION > 5
      filter->SetInputData(data2);
#else
      filter->SetInput(data2);
#endif
      filter->Update();
      data2 = filter->GetOutput();
    }
    UpdateFilteredImage(data1, data2);
  }
  OnLabelDataUpdated();
}

void LayerROI::Erode(int nTimes)
{
  SaveForUndo();
  if (m_layerMappedSurface)
  {
    ::LabelErode(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes);
  }
  else
  {
    vtkSmartPointer<vtkImageData> data1 = GetThresholdedMaskImage();
    vtkSmartPointer<vtkImageData> data2 = data1;
    for (int i = 0; i < nTimes; i++)
    {
      vtkSmartPointer<vtkImageDilateErode3D> filter = vtkSmartPointer<vtkImageDilateErode3D>::New();
      filter->SetKernelSize(2,2,2);
      filter->SetErodeValue(2);
      filter->SetDilateValue(0);
#if VTK_MAJOR_VERSION > 5
      filter->SetInputData(data2);
#else
      filter->SetInput(data2);
#endif
      filter->Update();
      data2 = filter->GetOutput();
    }
    UpdateFilteredImage(data1, data2);
  }
  OnLabelDataUpdated();
}

void LayerROI::Open(int nTimes)
{
  SaveForUndo();
  if (m_layerMappedSurface)
  {
    ::LabelErode(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes);
    ::LabelDilate(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes,
                  m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
  }
  else
  {
    vtkSmartPointer<vtkImageData> data1 = GetThresholdedMaskImage();
    vtkSmartPointer<vtkImageData> data2 = data1;
    for (int i = 0; i < nTimes; i++)
    {
      vtkSmartPointer<vtkImageOpenClose3D> filter = vtkSmartPointer<vtkImageOpenClose3D>::New();
      filter->SetKernelSize(2,2,2);
      filter->SetOpenValue(2);
      filter->SetCloseValue(0);
  #if VTK_MAJOR_VERSION > 5
      filter->SetInputData(data2);
  #else
      filter->SetInput(data2);
  #endif
      filter->Update();
      data2 = filter->GetOutput();
    }
    UpdateFilteredImage(data1, data2);
  }
  OnLabelDataUpdated();
}

void LayerROI::Close(int nTimes)
{
  SaveForUndo();
  if (m_layerMappedSurface)
  {
    ::LabelDilate(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes,
                  m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    ::LabelErode(m_label->GetRawLabel(), m_layerMappedSurface->GetSourceSurface()->GetMRIS(), nTimes);
  }
  else
  {
    vtkSmartPointer<vtkImageData> data1 = GetThresholdedMaskImage();
    vtkSmartPointer<vtkImageData> data2 = data1;
    for (int i = 0; i < nTimes; i++)
    {
      vtkSmartPointer<vtkImageOpenClose3D> filter = vtkSmartPointer<vtkImageOpenClose3D>::New();
      filter->SetKernelSize(2,2,2);
      filter->SetCloseValue(2);
      filter->SetOpenValue(0);
  #if VTK_MAJOR_VERSION > 5
      filter->SetInputData(data2);
  #else
      filter->SetInput(data2);
  #endif
      filter->Update();
      data2 = filter->GetOutput();
    }
    UpdateFilteredImage(data1, data2);
  }
  OnLabelDataUpdated();
}

void LayerROI::Resample()
{
  SaveForUndo();
  if (m_layerMappedSurface)
  {
    LABEL* old_label = m_label->GetRawLabel();
    ::LabelUnassign(old_label);
    LABEL* label = ::LabelSampleToSurface(m_layerMappedSurface->GetSourceSurface()->GetMRIS(), old_label,
                                          m_layerSource->GetSourceVolume()->GetMRI(),
                                          m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    if (label)
    {
      LabelCopy(label, old_label) ;
      LabelFree(&label);
    }
  }
  OnLabelDataUpdated();
}

bool LayerROI::HasUndo()
{
  return m_label->HasUndo();
}

bool LayerROI::HasRedo()
{
  return m_label->HasRedo();
}

void LayerROI::Undo()
{
  m_label->Undo();
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume() );
  SetModified();
  if (m_layerMappedSurface)
  {
    m_label->Initialize(m_layerSource->GetSourceVolume(), m_layerMappedSurface->GetSourceSurface(),
                        m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    m_layerMappedSurface->UpdateOverlay(true, true);
  }
  emit ActorUpdated();
  emit Modified();
}

void LayerROI::Redo()
{
  m_label->Redo();
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume() );
  SetModified();
  if (m_layerMappedSurface)
  {
    m_label->Initialize(m_layerSource->GetSourceVolume(), m_layerMappedSurface->GetSourceSurface(),
                        m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    m_layerMappedSurface->UpdateOverlay(true, true);
  }
  emit ActorUpdated();
  emit Modified();
}

void LayerROI::SaveForUndo(int nPlane, bool bAllFrame)
{
  Q_UNUSED(nPlane);
  Q_UNUSED(bAllFrame);
  m_label->SaveForUndo();
}

void LayerROI::Clear()
{
  m_label->Clear();
  m_label->UpdateRASImage( m_imageData, m_layerSource->GetSourceVolume() );
  SetModified();
  if (m_layerMappedSurface)
  {
    m_label->Initialize(m_layerSource->GetSourceVolume(), m_layerMappedSurface->GetSourceSurface(),
                        m_layerMappedSurface->IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    m_layerMappedSurface->UpdateOverlay(true, true);
  }
  emit ActorUpdated();
  emit Modified();
}

LABEL* LayerROI::GetRawLabel()
{
  return m_label->GetRawLabel();
}
