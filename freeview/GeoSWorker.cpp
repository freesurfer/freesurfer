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

GeoSWorker::GeoSWorker(QObject *parent) : QObject(parent)
{
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
}

void GeoSWorker::Compute(LayerMRI *mri, LayerMRI* seg, LayerMRI* seeds, int max_distance)
{
  m_mri = mri;
  m_seg = seg;
  m_seeds = seeds;
  if (max_distance > 0)
    m_nMaxDistance = max_distance;
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
  vtkSmartPointer<vtkExtractVOI> voi = vtkSmartPointer<vtkExtractVOI>::New();
  voi->SetInputConnection(cast->GetOutputPort());
  voi->SetVOI(bound);
  voi->Update();
  vtkSmartPointer<vtkImageCast> cast2 = vtkSmartPointer<vtkImageCast>::New();
#if VTK_MAJOR_VERSION > 5
  cast2->SetInputData(m_mri->GetImageData());
#else
  cast2->SetInput(m_mri->GetImageData());
#endif
  cast2->SetOutputScalarTypeToDouble();
  vtkSmartPointer<vtkExtractVOI> voi2 = vtkSmartPointer<vtkExtractVOI>::New();
  voi2->SetInputConnection(cast2->GetOutputPort());
  voi2->SetVOI(bound);
  voi2->Update();
  seeds = voi->GetOutput();
  vtkImageData* mri = voi2->GetOutput();
  GeodesicMatting geos;
  std::vector<unsigned char> label_list;
  label_list.push_back(1);
  label_list.push_back(2);
  double mri_range[2];
  vtkDoubleArray::SafeDownCast(mri->GetPointData()->GetScalars())->GetValueRange(mri_range);
  int dim_new[3];
  mri->GetDimensions(dim_new);
  vol_size = dim_new[0]*dim_new[1]*dim_new[2];
  unsigned char* seeds_out = new unsigned char[vol_size];
  bool bSuccess = geos.Compute(dim_new, (double*)mri->GetScalarPointer(), mri_range, (unsigned char*)seeds->GetScalarPointer(), label_list, seeds_out);
  if (bSuccess)
  {
//    m_seg->SaveForUndo();
    void* p = m_seg->GetImageData()->GetScalarPointer();
    int nDataType = m_seg->GetImageData()->GetScalarType();
    double fillValue = m_seg->GetFillValue();
    for (size_t n = 0; n < vol_size; n++)
    {
      if (seeds_out[n] > 0)
      {
        size_t i = (n%dim_new[0]), j = ((n/dim_new[0])%dim_new[1]), k = n/(dim_new[0]*dim_new[1]);
        i = (i+bound[0]) + (j+bound[2])*dim[0] + (k+bound[4])*dim[0]*dim[1];
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
  }
  else
    emit Failed();
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
