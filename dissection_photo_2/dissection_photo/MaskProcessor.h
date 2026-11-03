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
  // Shows components saved earlier (a SaveToNpy() file) on the loaded photo, as if
  // their boxes had been drawn: they are taken out of what new boxes can pick.
  bool LoadSaved(const QString& fn);
  // The colours GetMaskImage() gives a saved label file (empty when unreadable).
  static QImage SavedMaskImage(const QString& fn, const QList<QColor>& colors);

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
