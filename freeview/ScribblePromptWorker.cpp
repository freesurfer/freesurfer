#include "ScribblePromptWorker.h"
#include "TorchScriptModule.h"
#include "LayerMRI.h"
#include "vtkImageData.h"
#include "vtkImageCast.h"
#include "vtkImageShiftScale.h"
#include "vtkImageCast.h"
#include <QElapsedTimer>
#include <QDebug>
#include <QProcessEnvironment>
#include <QFile>

ScribblePromptWorker::ScribblePromptWorker(QObject *parent)
  : QObject(parent)
{
  m_module = new TorchScriptModule;
  connect(this, SIGNAL(ComputeTriggered()), SLOT(DoCompute()));
  connect(this, SIGNAL(ApplyTriggered()), SLOT(DoApply()));
  connect(this, SIGNAL(InitializationTriggered(QString)), SLOT(DoInitialization(QString)));
  QString fn = QProcessEnvironment::systemEnvironment().value( "FREESURFER_HOME" ) + "/traced_ScribblePrompt_UNet_nf192_res128.pt";
  if (QFile::exists(fn))
    Initialize(fn);
  else
    qDebug() << "Could not locate module file";
}

ScribblePromptWorker::~ScribblePromptWorker()
{
  m_module->deleteLater();
}

void ScribblePromptWorker::DoInitialization(const QString &fn)
{
  m_module->Load(fn);
}

void ScribblePromptWorker::Compute(LayerMRI *mri, LayerMRI* seg, LayerMRI* seeds, int nPlane, int nSlice, double fill_val)
{
  m_mri = mri;
  m_seg = seg;
  m_seeds = seeds;
  m_nInputPlane = nPlane;
  m_nInputSlice = nSlice;
  m_dFillValue = fill_val;
  emit ComputeTriggered();
}

void ScribblePromptWorker::Apply(LayerMRI *seg, LayerMRI *filled)
{
  m_seg = seg;
  m_filled = filled;
  emit ApplyTriggered();
}

void ScribblePromptWorker::DoCompute()
{
  QElapsedTimer timer;
  timer.start();
  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  cast->SetInputData(m_mri->GetSliceImageData(m_nInputPlane));
  cast->SetOutputScalarTypeToFloat();
  cast->Update();
  vtkImageData* img_mri = cast->GetOutput();
  int* dim = img_mri->GetDimensions();
  float* mri_ptr = (float*)img_mri->GetScalarPointer();
  unsigned char* seeds_ptr = (unsigned char*)m_seeds->GetSliceImageData(m_nInputPlane)->GetScalarPointer();
  bool bOverSize = (dim[0] > 128 || dim[1] > 128);
  int x_range[2] = {1000000,-1000000}, y_range[2] = {1000000,-1000000};
  int start_x = 0, start_y = 0;
  if (bOverSize)
  {
    for (int i = 0; i < dim[0]; i++)
    {
      for (int j = 0; j < dim[1]; j++)
      {
        if (seeds_ptr[j*dim[0]+i] > 0)
        {
          if (i < x_range[0])
            x_range[0] = i;
          if (i > x_range[1])
            x_range[1] = i;
          if (j < y_range[0])
            y_range[0] = j;
          if (j > y_range[1])
            y_range[1] = j;
        }
      }
    }

    if (x_range[1]-x_range[0] < 0 || y_range[1]-y_range[0] < 0)
    {
      qDebug() << "No foreground or background seeds selected";
      emit ComputeFinished(timer.elapsed()/1000.0);
      return;
    }
    else if (x_range[1]-x_range[0] > 128  || y_range[1]-y_range[0] > 128)
    {
      qDebug() << "FOV over 128 x 128";
      emit ComputeFinished(timer.elapsed()/1000.0);
      return;
    }
    start_x = qMax(0, (x_range[1]+x_range[0])/2-64);
    start_y = qMax(0, (y_range[1]+y_range[0])/2-64);
  }

  QVector<float*> inputs;
  for (int i = 0; i < 5; i++)
  {
    float* ptr = new float[128*128];
    memset(ptr, 0, sizeof(float)*128*128);
    inputs << ptr;
  }

  double value_r[2];
  img_mri->GetScalarRange(value_r);
  for (int i = 0; i < 128; i++)
  {
    for (int j = 0; j < 128; j++)
    {
      int x = start_x + i, y = start_y + j;
      if (x >= dim[0] || y >= dim[1])
        continue;

      inputs[0][j*128+i] = mri_ptr[y*dim[0]+x]/value_r[1];
      if (seeds_ptr[y*dim[0]+x] == 1)       // foreground
        inputs[2][j*128+i] = 1;
      else if (seeds_ptr[y*dim[0]+x] == 2)    // background
        inputs[3][j*128+i] = 1;
      else if (seeds_ptr[y*dim[0]+x] == 3)   // box
        inputs[1][j*128+i] = 1;
    }
  }

  float* output = new float[128*128];
  m_module->Run(inputs, output);
  void* p = m_seg->GetImageData()->GetScalarPointer();
  int nDataType = m_seg->GetImageData()->GetScalarType();
  double fillValue = m_seg->GetFillValue();
  dim = m_seg->GetImageData()->GetDimensions();
  int x, y;
  for (int i = 0; i < 128; i++)
  {
    for (int j = 0; j < 128; j++)
    {
      if (output[j*128+i] <= 0)
        continue;

      int x = start_x + i, y = start_y + j;
      int n = 0;
      switch (m_nInputPlane)
      {
      case 0:
        if (y >= dim[2] || x >= dim[1])
          continue;
        n = y*dim[1]*dim[0] + x*dim[0] + m_nInputSlice;
        break;
      case 1:
        if (y >= dim[2] || x >= dim[0])
          continue;
        n = y*dim[1]*dim[0] + m_nInputSlice*dim[0] + x;
        break;
      case 2:
        if (y >= dim[1] || x >= dim[0])
          continue;
        n = m_nInputSlice*dim[1]*dim[0] + y*dim[0] + x;
        break;
      }

      switch (nDataType)
      {
      case VTK_INT:
        ((int*)p)[n] = (int)fillValue;
        break;
      case VTK_UNSIGNED_CHAR:
        ((unsigned char*)p)[n] = (unsigned char)fillValue;
        break;
      case VTK_FLOAT:
        ((float*)p)[n] = (float)fillValue;
        break;
      case VTK_DOUBLE:
        ((double*)p)[n] = (double)fillValue;
        break;
      }
    }
  }

  for (int i = 0; i < inputs.size(); i++)
    delete[] inputs[i];
  delete[] output;

  m_seg->SetModified();
  emit ComputeFinished(timer.elapsed()/1000.0);
}

void ScribblePromptWorker::DoApply()
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
