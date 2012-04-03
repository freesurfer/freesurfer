#include "LayerMRIWorkerThread.h"
#include "LayerMRI.h"
#include "vtkImageData.h"
#include <QMutexLocker>

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

  IntList vals;
  for (int i = 0; i < dim[0]; i++)
  {
    for (int j = 0; j < dim[1]; j++)
    {
      for (int k = 0; k < dim[2]; k++)
      {
        int val = (int)image->GetScalarComponentAsDouble(i, j, k, 0);
        if (val != 0 && !vals.contains(val))
          vals << val;
      }
      {
        QMutexLocker locker(&mutex);
        if (m_bAbort)
          return;
      }
    }
  }
  emit AvailableLabels(vals);
}
