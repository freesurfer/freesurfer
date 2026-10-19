#include "MaskProcessor.h"
#include <QDebug>
#include <QFile>
#include "cnpy.h"

MaskProcessor::MaskProcessor() : m_data(NULL), m_dataInBuffer(NULL), m_dataOutBuffer(NULL),
  m_nWidth(0), m_nHeight(0)
{
}

MaskProcessor::~MaskProcessor()
{
  ClearData();
}

void MaskProcessor::ClearData()
{
  delete[] m_data;
  m_data = NULL;
  delete[] m_dataInBuffer;
  m_dataInBuffer = NULL;
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
  if (!m_data || !m_dataOutBuffer)
    return false;

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
  for (int i = 0; i + 1 < pts.size(); i+=2)
  {
    ProcessSelection(pts[i], pts[i+1], i/2+1);
  }
}

QImage MaskProcessor::GetMaskImage(const QList<QColor>& colors)
{
  if (!m_dataOutBuffer || colors.isEmpty())
    return QImage();

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

bool MaskProcessor::LoadSaved(const QString& fn)
{
  if (!m_data || !m_dataOutBuffer || !QFile::exists(fn))
    return false;
  cnpy::NpyArray ar;
  try { ar = cnpy::npy_load(qPrintable(fn)); }
  catch (const std::exception&) { return false; }
  if (ar.shape.size() != 2 || (int)ar.shape[0] != m_nHeight || (int)ar.shape[1] != m_nWidth || ar.word_size != 1)
    return false;
  const unsigned char* p = ar.data<unsigned char>();
  for (int n = 0; n < m_nWidth*m_nHeight; n++)
  {
    if (p[n] > 0)
    {
      m_dataOutBuffer[n] = p[n];
      m_data[n] = 0;
    }
  }
  return true;
}

QImage MaskProcessor::SavedMaskImage(const QString& fn, const QList<QColor>& colors)
{
  if (!QFile::exists(fn) || colors.isEmpty())
    return QImage();
  cnpy::NpyArray ar;
  try { ar = cnpy::npy_load(qPrintable(fn)); }
  catch (const std::exception&) { return QImage(); }
  if (ar.shape.size() != 2 || ar.word_size != 1)
    return QImage();
  const int h = ar.shape[0], w = ar.shape[1];
  const unsigned char* p = ar.data<unsigned char>();
  QImage image(w, h, QImage::Format_ARGB32);
  image.fill(QColor(0,0,0,0));
  for (int j = 0; j < h; j++)
  {
    QRgb* q = (QRgb*)image.scanLine(j);
    for (int i = 0; i < w; i++)
      if (p[j*w+i] > 0)
        q[i] = colors[(p[j*w+i]-1)%colors.size()].rgb();
  }
  return image;
}

bool MaskProcessor::SaveToNpy(const QString& fn)
{
  std::vector<size_t> shape;
  shape.push_back(m_nHeight);
  shape.push_back(m_nWidth);
  cnpy::npy_save<unsigned char>(qPrintable(fn), m_dataOutBuffer, shape, "w");
  return true;
}
