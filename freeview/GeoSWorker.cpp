#include "GeoSWorker.h"
#include "geos/GeodesicMatting.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "LayerMRI.h"
#include "LayerROI.h"
#include "vtkImageWeightedSum.h"
#include "vtkImageCast.h"
#include "vtkImageDilateErode3D.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include <QDebug>
#include <QFile>

GeoSWorker::GeoSWorker(QObject *parent) : QObject(parent)
{
  connect(this, SIGNAL(ComputeTriggered()), SLOT(DoCompute()));
  moveToThread(&m_thread);
  m_thread.start();
}

GeoSWorker::~GeoSWorker()
{
  m_thread.quit();
  m_thread.wait();
}

void GeoSWorker::Compute(LayerMRI *mri, LayerMRI* seg, LayerMRI* seeds)
{
  m_mri = mri;
  m_seg = seg;
  m_seeds = seeds;
  emit ComputeTriggered();
}

void GeoSWorker::DoCompute()
{
  int* dim = m_seg->GetImageData()->GetDimensions();
  size_t vol_size = dim[0]*dim[1]*dim[2];
//  vtkSmartPointer<vtkImageWeightedSum> sum = vtkSmartPointer<vtkImageWeightedSum>::New();
//  sum->AddInput(m_roiInterior->GetImageData());
//  sum->AddInput(m_roiExterior->GetImageData());
//  sum->SetWeight(0, 1);
//  sum->SetWeight(1, 1);

//  float* ptr1 = (float*)m_roiInterior->GetImageData()->GetScalarPointer();
//  float* ptr2 = (float*)m_roiExterior->GetImageData()->GetScalarPointer();
//  for (int i = 0; i < vol_size; i++)
//    ptr1[i] += ptr2[i];

  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  cast->SetInput(m_seg->GetImageData());
  cast->SetOutputScalarTypeToUnsignedChar();
  cast->Update();
  vtkImageData* seeds = cast->GetOutput();
//  QFile file("/autofs/space/voxel_001/users/rpwang/src/GeodesicMatting/foo.img");
//  file.open(QFile::WriteOnly);
//  file.write((const char*)seeds->GetScalarPointer(), sizeof(char)*vol_size);
//  file.close();

  // find the VOI in seeds
  //  int bound[6] = {0, dim[0], 0, dim[1], 0, dim[2]};
  //  unsigned char* ptr = (unsigned char*)seeds->GetScalarPointer();
  //  for (int i = 0; )
  //  vtkSmartPointer<vtkExtractVOI> voi = vtkSmartPointer<vtkExtractVOI>::New();
  //  voi->SetInput(seeds);
  //  voi->SetVOI(bound);
  vtkSmartPointer<vtkImageCast> cast2 = vtkSmartPointer<vtkImageCast>::New();
  cast2->SetInput(m_mri->GetImageData());
  cast2->SetOutputScalarTypeToDouble();
  cast2->Update();
  vtkImageData* mri = cast2->GetOutput();
  GeodesicMatting geos;
  std::vector<unsigned char> label_list;
  label_list.push_back(1);
  label_list.push_back(2);
  double mri_range[2];
  vtkDoubleArray::SafeDownCast(mri->GetPointData()->GetScalars())->GetValueRange(mri_range);
  unsigned char* seeds_out = new unsigned char[dim[0]*dim[1]*dim[2]];
  bool bSuccess = geos.ComputeWithSorting(dim, (double*)mri->GetScalarPointer(), mri_range, (unsigned char*)seeds->GetScalarPointer(), label_list, seeds_out);
  if (bSuccess)
  {
    void* p = m_seg->GetImageData()->GetScalarPointer();
    int nDataType = m_seg->GetImageData()->GetScalarType();
    double fillValue = m_seg->GetFillValue();
    for (size_t i = 0; i < vol_size; i++)
    {
      if (seeds_out > 0)
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
    m_seg->Modified();
  }
  emit Finished(bSuccess);
}
