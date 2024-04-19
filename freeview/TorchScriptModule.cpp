#include "TorchScriptModule.h"
#undef slots
#undef Byte
#include "torch/script.h"
#define slots Q_SLOTS
#include <QDebug>

TorchScriptModule::TorchScriptModule(QObject *parent)
    : QObject(parent)
{
  m_module = new torch::jit::script::Module;
}

TorchScriptModule::~TorchScriptModule()
{
  delete ((torch::jit::script::Module*)m_module);
}

void TorchScriptModule::Load(const QString &fn)
{
  *((torch::jit::script::Module*)m_module) = torch::jit::load(qPrintable(fn));
}

void FillInTensorFromBuffer(torch::Tensor &t, int n, float *buf_in)
{
  auto ptr = t.accessor<float,4>();
  for (int i = 0; i < 128; i++)
  {
    for (int j = 0; j < 128; j++)
    {
      ptr[0][n][j][i] = buf_in[j*128+i];
    }
  }
}

void FillInBufferFromTensor(float *buffer, torch::Tensor &t_in, int n = 0)
{
  auto tptr = t_in.accessor<float,4>();
  for (int i = 0; i < 128; i++)
  {
    for (int j = 0; j < 128; j++)
    {
      buffer[j*128+i] = tptr[0][n][j][i];
    }
  }
}

void TorchScriptModule::Run(QVector<float*> in_ptr, float *output)
{
  try {
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor t = torch::zeros({1, in_ptr.size(), 128, 128});
    for (int i = 0; i < in_ptr.size(); i++)
    {
      if (in_ptr[i])
        FillInTensorFromBuffer(t, i, in_ptr[i]);
    }
    inputs.push_back(t);
    // Execute the model and turn its output into a tensor.
    torch::Tensor t_out = ((torch::jit::script::Module*)m_module)->forward(inputs).toTensor();
    FillInBufferFromTensor(output, t_out);
  }
  catch (const c10::Error& e) {
    qDebug() << "error running the module\n";
  }
  emit Finished();
}

