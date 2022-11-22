#ifndef MASKPROCESSOR_H
#define MASKPROCESSOR_H

#include <QString>
#include <QPoint>
#include <QList>
#include <QColor>
#include <QImage>

#ifndef PY_DATA_TYPE
#define PY_DATA_TYPE uint16_t
#endif

class MaskProcessor
{
public:
  MaskProcessor();  
  ~MaskProcessor();

  bool Load(const QString& mask_fn);
  void LoadSelections(const QList<QPoint>& pts);
  bool ProcessSelection(const QPoint& pt1, const QPoint& pt2, int nFillValue);
  bool SaveToNpy(const QString& fn);

  void ClearBuffer();

  QImage GetMaskImage(const QList<QColor>& colors);

private:
  void ClearData();

  PY_DATA_TYPE* m_data;
  PY_DATA_TYPE* m_dataInBuffer;
  unsigned char* m_dataOutBuffer;
  int m_nWidth;
  int m_nHeight;
};

#endif // MASKPROCESSOR_H
