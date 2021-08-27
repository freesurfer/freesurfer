#include "LayerMRIWorkerThread.h"
#include "LayerMRI.h"
#include "vtkImageData.h"
#include <QMutexLocker>
#include "MyVTKUtils.h"
#include <QDateTime>

LayerMRIWorkerThread::LayerMRIWorkerThread(LayerMRI *mri) :
  QThread(mri), m_bAbort(false)
{
}

void LayerMRIWorkerThread::Abort()
{
  QMutexLocker locker(&mutex);
  m_bAbort = true;
}

void LayerMRIWorkerThread::run()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(parent());
  vtkImageData* image = mri->GetImageData();
  int* dim = image->GetDimensions();
  double* origin = image->GetOrigin();
  double* vs = image->GetSpacing();

  IntList vals;
  QMap<int, QList<double> > centers;
  QMap<int, int> counts;
  char* ptr = (char*)image->GetScalarPointer();
  int scalar_type = image->GetScalarType();
  int n_frames = image->GetNumberOfScalarComponents();
  for (int i = 0; i < dim[0]; i++)
  {
    for (int j = 0; j < dim[1]; j++)
    {
      for (int k = 0; k < dim[2]; k++)
      {
        int val = (int)MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, j, k, 0, scalar_type);
        if (val != 0)
        {
          if (!vals.contains(val))
            vals << val;
          if (!centers.contains(val))
          {
            QList<double> center;
            center << i*vs[0] + origin[0] << j*vs[1] + origin[1] << k*vs[2] + origin[2];
            centers[val] = center;
            counts[val] = 1;
          }
          else
          {
            centers[val][0] += i*vs[0] + origin[0];
            centers[val][1] += j*vs[1] + origin[1];
            centers[val][2] += k*vs[2] + origin[2];
            counts[val] ++;
          }
        }
      }
      {
        QMutexLocker locker(&mutex);
        if (m_bAbort)
          return;
      }
    }
  }
  QList<int> keys = centers.keys();
  foreach(int val, keys)
  {
    centers[val][0] /= counts[val];
    centers[val][1] /= counts[val];
    centers[val][2] /= counts[val];
  }

  QMutexLocker locker(&mutex);
  mri->m_nAvailableLabels = vals;
  mri->m_listLabelCenters = centers;
  mri->m_labelVoxelCounts = counts;
  mri->setProperty("stats_last_updated", QDateTime::currentMSecsSinceEpoch());

  emit LabelInformationReady();
}
