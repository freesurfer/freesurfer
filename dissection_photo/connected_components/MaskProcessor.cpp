#include "MaskProcessor.h"
#include <QDebug>
#include <QFile>
#include "cnpy.h"
#include <QList>

MaskProcessor::MaskProcessor() : m_data(NULL), m_dataInBuffer(NULL), m_dataOutBuffer(NULL)
{

}

MaskProcessor::~MaskProcessor()
{
  ClearData();
}

void MaskProcessor::ClearData()
{
  if (m_data)
    delete[] m_data;
  m_data = NULL;

  if (m_dataInBuffer)
    delete[] m_dataInBuffer;
  m_dataInBuffer = NULL;

  if (m_dataOutBuffer)
    delete[] m_dataOutBuffer;
  m_dataOutBuffer = NULL;
}

bool MaskProcessor::Load(const QString& npy_in)
{
  if (!QFile::exists(npy_in))
  {
    qDebug() << "File does not exist:" << npy_in;
    return false;
  }
  cnpy::NpyArray ar = cnpy::npz_load(qPrintable(npy_in), "cc");
  if (ar.shape.size() != 2)
  {
    qDebug() << "Could not load numpy file " << npy_in;
    return false;
  }

  ClearData();
  PY_DATA_TYPE* ptr = ar.data<PY_DATA_TYPE>();
  int nsize = ar.shape[1] * ar.shape[0];
  m_nWidth = ar.shape[1];
  m_nHeight = ar.shape[0];

  m_dataInBuffer = new PY_DATA_TYPE[nsize];
  memcpy(m_dataInBuffer, ptr, nsize*sizeof(PY_DATA_TYPE));

  m_data = new PY_DATA_TYPE[nsize];
  memcpy(m_data, ptr, nsize*sizeof(PY_DATA_TYPE));

  m_dataOutBuffer = new unsigned char[nsize];
  memset(m_dataOutBuffer, 0, nsize);
  return true;
}

bool MaskProcessor::ProcessSelection(const QPoint& pt1, const QPoint& pt2, int nFillValue)
{
  QList<int> list;
  for (int i = pt1.x(); i <= pt2.x(); i++)
  {
    for (int j = pt1.y(); j <= pt2.y(); j++)
    {
      if (i >= 0 && i < m_nWidth && j >= 0 && j < m_nHeight)
      {
        int nVal = m_data[j*m_nWidth+i];
        if (nVal > 0 && !list.contains(nVal))
          list << nVal;
      }
    }
  }
  if (list.isEmpty())
    return false;

  for (int i = 0; i < m_nWidth; i++)
  {
    for (int j = 0; j < m_nHeight; j++)
    {
      int nIndex = j*m_nWidth+i;
      if (list.contains(m_data[nIndex]))
      {
        m_dataOutBuffer[nIndex] = nFillValue;
        m_data[nIndex] = 0;
      }
    }
  }
  return true;
}

void MaskProcessor::LoadSelections(const QList<QPoint>& pts)
{
  for (int i = 0; i < pts.size(); i+=2)
  {
    ProcessSelection(pts[i], pts[i+1], i/2+1);
  }
}

QImage MaskProcessor::GetMaskImage(const QList<QColor>& colors)
{
  QImage image(m_nWidth, m_nHeight, QImage::Format_ARGB32);
  image.fill(QColor(0,0,0,0));
  for (int j = 0; j < image.height(); j++)
  {
    QRgb* p = (QRgb*)image.scanLine(j);
    for (int i = 0; i < image.width(); i++)
    {
      int n = j*image.width()+i;
      if (m_dataOutBuffer[n] > 0)
      {
        p[i] = colors[(m_dataOutBuffer[n]-1)%colors.size()].rgb();
      }
    }
  }
  return image;
}

void MaskProcessor::ClearBuffer()
{
  if (!m_data || !m_dataOutBuffer)
    return;

  memcpy(m_data, m_dataInBuffer, m_nWidth*m_nHeight*sizeof(PY_DATA_TYPE));
  memset(m_dataOutBuffer, 0, m_nWidth*m_nHeight);
}

bool MaskProcessor::SaveToNpy(const QString& fn)
{
  std::vector<size_t> shape;
  shape.push_back(m_nHeight);
  shape.push_back(m_nWidth);
  cnpy::npy_save<unsigned char>(qPrintable(fn), m_dataOutBuffer, shape, "w");
  return true;
}
