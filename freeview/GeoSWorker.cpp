#include "GeoSWorker.h"
#include "geos/GeodesicMatting.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "LayerMRI.h"
#include "LayerROI.h"
#include "vtkImageCast.h"
#include "vtkImageDilateErode3D.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkExtractVOI.h"
#include <QDebug>
#include <QFile>
#include "vtkImageGaussianSmooth.h"
#include "vtkImageExtractComponents.h"
#include "vtkImageThreshold.h"
#include <QElapsedTimer>
#include "LayerPropertyMRI.h"

GeoSWorker::GeoSWorker(QObject *parent) : QObject(parent)
{
  m_geos = new GeodesicMatting;
  connect(m_geos, SIGNAL(Progress(double)), this, SIGNAL(Progress(double)));
  m_nMaxDistance = 10;
  connect(this, SIGNAL(ComputeTriggered()), SLOT(DoCompute()));
  connect(this, SIGNAL(ApplyTriggered()), SLOT(DoApply()));
  moveToThread(&m_thread);
  m_thread.start();
}

GeoSWorker::~GeoSWorker()
{
  m_thread.quit();
  m_thread.wait();
  m_geos->deleteLater();
}

void GeoSWorker::Compute(LayerMRI *mri, LayerMRI* seg, LayerMRI* seeds, int max_distance, double smoothing, LayerMRI* mask, double fill_val, int max_foreground_dist)
{
  m_mri = mri;
  m_seg = seg;
  m_seeds = seeds;
  m_dSmoothing = smoothing;
  m_mask = mask;
  m_nMaxForegroundDistance = max_foreground_dist;
  if (max_distance > 0)
    m_nMaxDistance = max_distance;
  m_dFillValue = fill_val;
  emit ComputeTriggered();
}

void GeoSWorker::Apply(LayerMRI *seg, LayerMRI *filled)
{
  m_seg = seg;
  m_filled = filled;
  emit ApplyTriggered();
}

void GeoSWorker::DoCompute()
{
  int* dim = m_seg->GetImageData()->GetDimensions();
  size_t vol_size = dim[0]*dim[1]*dim[2];

  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
#if VTK_MAJOR_VERSION > 5
  cast->SetInputData(m_seeds->GetImageData());
#else
  cast->SetInput(m_seeds->GetImageData());
#endif
  cast->SetOutputScalarTypeToUnsignedChar();
  cast->Update();
  vtkImageData* seeds = cast->GetOutput();
  // find the VOI in seeds
  int bound[6] = {dim[0], 0, dim[1], 0, dim[2], 0};
  unsigned char* ptr = (unsigned char*)seeds->GetScalarPointer();
  for (int i = 0; i < dim[0]; i++)
  {
    for (int j = 0; j < dim[1]; j++)
    {
      for (int k = 0; k < dim[2]; k++)
      {
        int nVal = ptr[k*dim[0]*dim[1]+j*dim[0]+i];
        if (nVal > 1)
        {
          if (i < bound[0])
            bound[0] = i;
          if (i > bound[1])
            bound[1] = i;
          if (j < bound[2])
            bound[2] = j;
          if (j > bound[3])
            bound[3] = j;
          if (k < bound[4])
            bound[4] = k;
          if (k > bound[5])
            bound[5] = k;
        }
        else if (nVal == 1)
        {
          if (i-m_nMaxDistance < bound[0])
            bound[0] = i-m_nMaxDistance;
          if (i+m_nMaxDistance > bound[1])
            bound[1] = i+m_nMaxDistance;
          if (j-m_nMaxDistance < bound[2])
            bound[2] = j-m_nMaxDistance;
          if (j+m_nMaxDistance > bound[3])
            bound[3] = j+m_nMaxDistance;
          if (k-m_nMaxDistance < bound[4])
            bound[4] = k-m_nMaxDistance;
          if (k+m_nMaxDistance > bound[5])
            bound[5] = k+m_nMaxDistance;
        }
      }
    }
  }
  bound[0] = qMax(0, bound[0]);
  bound[1] = qMin(dim[0]-1, bound[1]);
  bound[2] = qMax(0, bound[2]);
  bound[3] = qMin(dim[1]-1, bound[3]);
  bound[4] = qMax(0, bound[4]);
  bound[5] = qMin(dim[2]-1, bound[5]);
  if (bound[0] > bound[1])
    memset(bound, 0, sizeof(int)*6);
  vtkSmartPointer<vtkExtractVOI> voi = vtkSmartPointer<vtkExtractVOI>::New();
  voi->SetInputConnection(cast->GetOutputPort());
  voi->SetVOI(bound);
  voi->Update();
  vtkSmartPointer<vtkImageData> src_image = m_mri->GetImageData();
  vtkSmartPointer<vtkImageData> grayscale = src_image;
  // deal with RGB images
  if (m_mri->GetDataType() == MRI_RGB || (m_mri->GetNumberOfFrames() == 3 && m_mri->GetProperty()->GetDisplayRGB()))
  {
    grayscale = vtkSmartPointer<vtkImageData>::New();
    grayscale->SetSpacing(src_image->GetSpacing());
    grayscale->SetOrigin(src_image->GetOrigin());
    grayscale->SetDimensions(src_image->GetDimensions());
#if VTK_MAJOR_VERSION > 5
    grayscale->AllocateScalars(VTK_DOUBLE, 1);
#else
    grayscale->SetScalarTypeToDouble();
    grayscale->AllocateScalars();
#endif
    if (m_mri->GetDataType() == MRI_RGB)
    {
      unsigned char* in_ptr = (unsigned char*)src_image->GetScalarPointer();
      double* out_ptr = (double*)grayscale->GetScalarPointer();
      int* dim = src_image->GetDimensions();
      qlonglong nSize = ((qlonglong)dim[0])*dim[1]*dim[2];
      for (qlonglong i = 0; i < nSize; i++)
      {
        unsigned char* ptr = in_ptr + i*4;
        out_ptr[i] = ((double)ptr[0]) + ptr[1] + ptr[2];
      }
    }
    else
    {
      void* in_ptr = src_image->GetScalarPointer();
      double* out_ptr = (double*)grayscale->GetScalarPointer();
      int* dim = src_image->GetDimensions();
      qlonglong nSize = ((qlonglong)dim[0])*dim[1]*dim[2];
      switch (src_image->GetScalarType())
      {
      case VTK_UNSIGNED_CHAR:
        for (qlonglong i = 0; i < nSize; i++)
        {
          unsigned char* ptr = (unsigned char*)in_ptr + i*3;
          out_ptr[i] = ((double)ptr[0]) + ptr[1] + ptr[2];
        }
        break;
      case VTK_INT:
        for (qlonglong i = 0; i < nSize; i++)
        {
          int* ptr = (int*)in_ptr + i*3;
          out_ptr[i] = ((double)ptr[0]) + ptr[1] + ptr[2];
        }
        break;
      case VTK_SHORT:
        for (qlonglong i = 0; i < nSize; i++)
        {
          short* ptr = (short*)in_ptr + i*3;
          out_ptr[i] = ((double)ptr[0]) + ptr[1] + ptr[2];
        }
        break;
      case VTK_FLOAT:
        for (qlonglong i = 0; i < nSize; i++)
        {
          float* ptr = (float*)in_ptr + i*3;
          out_ptr[i] = ((double)ptr[0]) + ptr[1] + ptr[2];
        }
        break;
      case VTK_DOUBLE:
        for (qlonglong i = 0; i < nSize; i++)
        {
          double* ptr = (double*)in_ptr + i*3;
          out_ptr[i] = ptr[0] + ptr[1] + ptr[2];
        }
        break;
      }
    }
  }
  vtkSmartPointer<vtkImageCast> cast2 = vtkSmartPointer<vtkImageCast>::New();
#if VTK_MAJOR_VERSION > 5
  cast2->SetInputData(grayscale);
#else
  cast2->SetInput(grayscale);
#endif
  cast2->SetOutputScalarTypeToDouble();
  vtkSmartPointer<vtkExtractVOI> voi2 = vtkSmartPointer<vtkExtractVOI>::New();
  voi2->SetInputConnection(cast2->GetOutputPort());
  voi2->SetVOI(bound);
  voi2->Update();
  seeds = voi->GetOutput();
  vtkSmartPointer<vtkImageData> mri = voi2->GetOutput();
  if (m_dSmoothing > 0)
  {
    vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetDimensionality(3);
    smooth->SetRadiusFactor(3);
    smooth->SetStandardDeviation(m_dSmoothing);
    smooth->SetInputConnection(voi2->GetOutputPort());
    smooth->Update();
    mri = smooth->GetOutput();
  }

  std::vector<unsigned char> label_list;
  label_list.push_back(1);
  label_list.push_back(2);
  double mri_range[2];
  vtkDoubleArray::SafeDownCast(mri->GetPointData()->GetScalars())->GetValueRange(mri_range);
  int dim_new[3];
  mri->GetDimensions(dim_new);
  vol_size = dim_new[0]*dim_new[1]*dim_new[2];
  unsigned char* seeds_out = new unsigned char[vol_size];
  QElapsedTimer timer;
  timer.start();
  unsigned char* seed_ptr = (unsigned char*)seeds->GetScalarPointer();

  float* mask_ptr = NULL;
  if (m_mask)
  {
    cast = vtkSmartPointer<vtkImageCast>::New();
#if VTK_MAJOR_VERSION > 5
    cast->SetInputData(m_mask->GetImageData());
#else
    cast->SetInput(m_mask->GetImageData());
#endif
    cast->SetOutputScalarTypeToFloat();
    cast->Update();
    vtkImageData* mask_image = cast->GetOutput();
    mask_ptr = (float*)mask_image->GetScalarPointer();
    for (int i = 0; i < dim_new[0]; i++)
    {
      for (int j = 0; j < dim_new[1]; j++)
      {
        for (int k = 0; k < dim_new[2]; k++)
        {
          size_t n_voi = k*dim_new[1]*dim_new[0] + j*dim_new[0] + i;
          size_t n = (k+bound[4])*dim[1]*dim[0] + (j+bound[2])*dim[0] + i+bound[0];
          if (mask_ptr[n] != 0)
          {
            if (mask_ptr[n] == m_dFillValue)
              seed_ptr[n_voi] = 1;
            else
              seed_ptr[n_voi] = 2;
          }
        }
      }
    }
  }

  unsigned char* fg_mask_ptr = NULL;
  vtkSmartPointer<vtkImageData> fg_mask_image;
  if (m_nMaxForegroundDistance > 0)
  {
    vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
    threshold->ThresholdBetween(1, 1);
    threshold->SetInValue(1);
    threshold->SetOutValue(0);
    threshold->ReplaceInOn();
    threshold->ReplaceOutOn();
#if VTK_MAJOR_VERSION > 5
    threshold->SetInputData(seeds);
#else
    threshold->SetInput(seeds);
#endif
    vtkSmartPointer<vtkImageDilateErode3D> dilate = vtkSmartPointer<vtkImageDilateErode3D>::New();
    dilate->SetInputConnection(threshold->GetOutputPort());
    dilate->SetKernelSize(m_nMaxForegroundDistance*2, m_nMaxForegroundDistance*2, m_nMaxForegroundDistance*2);
    dilate->SetDilateValue(1);
    dilate->SetErodeValue(0);
    dilate->Update();
    fg_mask_image = dilate->GetOutput();
    fg_mask_ptr = (unsigned char*)fg_mask_image->GetScalarPointer();
  }

  double scale[3] = {1,1,1};
  bool bSuccess = m_geos->ComputeWithBinning(dim_new, scale, (double*)mri->GetScalarPointer(), mri_range, seed_ptr, label_list, seeds_out);
  if (bSuccess)
  {
    void* p = m_seg->GetImageData()->GetScalarPointer();
    int nDataType = m_seg->GetImageData()->GetScalarType();
    double fillValue = m_seg->GetFillValue();
    for (size_t n = 0; n < vol_size; n++)
    {
      if (seeds_out[n] > 0)
      {
        size_t i = (n%dim_new[0]), j = ((n/dim_new[0])%dim_new[1]), k = n/(dim_new[0]*dim_new[1]);
        i = (i+bound[0]) + (j+bound[2])*dim[0] + (k+bound[4])*dim[0]*dim[1];
        if ((!mask_ptr || mask_ptr[i] == 0) && (!fg_mask_ptr || fg_mask_ptr[n] > 0))
        {
          switch (nDataType)
          {
          case VTK_INT:
            ((int*)p)[i] = (int)fillValue;
            break;
          case VTK_UNSIGNED_CHAR:
            ((unsigned char*)p)[i] = (unsigned char)fillValue;
            break;
          case VTK_FLOAT:
            ((float*)p)[i] = (float)fillValue;
            break;
          case VTK_DOUBLE:
            ((double*)p)[i] = (double)fillValue;
            break;
          }
        }
      }
    }
    m_seg->SetModified();
    emit ComputeFinished(timer.elapsed()/1000.0);
  }
  else
    emit ComputeFinished(-1);
}

void GeoSWorker::DoApply()
{
  m_seg->SaveForUndo();
  void* p = m_seg->GetImageData()->GetScalarPointer();
  int nDataType = m_seg->GetImageData()->GetScalarType();
  double fillValue = m_seg->GetFillValue();
  int* dim = m_filled->GetImageData()->GetDimensions();
  unsigned char* p_filled = (unsigned char*)m_filled->GetImageData()->GetScalarPointer();
  size_t vol_size = dim[0]*dim[1]*dim[2];
  for (size_t i = 0; i < vol_size; i++)
  {
    if (p_filled[i])
    {
      switch (nDataType)
      {
      case VTK_INT:
        ((int*)p)[i] = (int)fillValue;
        break;
      case VTK_UNSIGNED_CHAR:
        ((unsigned char*)p)[i] = (unsigned char)fillValue;
        break;
      case VTK_FLOAT:
        ((float*)p)[i] = (float)fillValue;
        break;
      case VTK_DOUBLE:
        ((double*)p)[i] = (double)fillValue;
        break;
      }
    }
  }
  m_seg->SetModified();
  emit ApplyFinished();
}

void GeoSWorker::Abort()
{
  m_geos->Abort();
}

QString GeoSWorker::GetErrorMessage()
{
  return m_geos?m_geos->GetErrorMessage():"";
}
