#include "VolumeFilterOptimal.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "LayerROI.h"
#include "ProgressCallback.h"
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include <QList>
#include <vector>
#include <QDebug>

VolumeFilterOptimal::VolumeFilterOptimal(QList<LayerMRI*> inputs, QList<LayerROI*> input_labels, LayerMRI* output, QObject* parent) :
  VolumeFilter(inputs[0], output, parent),
  m_inputMRIs(inputs),
  m_inputROIs(input_labels)
{
}

bool VolumeFilterOptimal::Execute()
{
  if (m_inputROIs.size() < 2)
    return false;

  ::SetProgressCallback(ProgressCallback, 10, 50);
  std::vector<void*> images;
  vtkSmartPointer<vtkImageData> image;
  foreach (LayerMRI* mri, m_inputMRIs)
  {
    vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
    cast->SetInputData(mri->GetImageData());
    cast->SetOutputScalarTypeToFloat();
    cast->Update();
    image = cast->GetOutput();
    images.push_back((void*)image->GetScalarPointer());
  }

  QList<int> list[2];
  for (int n = 0; n < 2; n++)
  {
    vtkImageData* roi_data = m_inputROIs[n]->GetImageData();
    int* dim = roi_data->GetDimensions();
    float* ptr = (float*)roi_data->GetScalarPointer();
    for (int i = 0; i < dim[0]; i++)
    {
      for (int j = 0; j < dim[1]; j++)
      {
        for (int k = 0; k < dim[2]; k++)
        {
          int ndx = k*dim[0]*dim[1] + j*dim[0] + i;
          if (ptr[ndx] >= 0)
            list[n] << ndx;
        }
      }
    }
  }

  int* vox0 = new int[list[0].size()];
  int* vox1 = new int[list[1].size()];
  for (int i = 0; i < list[0].size(); i++)
    vox0[i] = list[0][i];
  for (int i = 0; i < list[1].size(); i++)
    vox1[i] = list[1][i];

  ::SetProgressCallback(ProgressCallback, 50, 60);
  float* out_ptr = (float*)m_volumeOutput->GetImageData()->GetScalarPointer();
  int* dim = m_volumeOutput->GetImageData()->GetDimensions();
  bool ret = MyUtils::CalculateOptimalVolume(vox0, list[0].size(), vox1, list[1].size(), images, out_ptr, dim[0]*dim[1]*dim[2]);
  ::SetProgressCallback(ProgressCallback, 90, 100);
  delete[] vox0;
  delete[] vox1;
  return ret;
}
