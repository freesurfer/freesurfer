#ifndef WIDGETIMAGEVIEW_H
#define WIDGETIMAGEVIEW_H

#include <QWidget>
#include <QImage>

#ifndef RECT_REGION
#define RECT_REGION QPair<QPoint,QPoint>
#endif

class WidgetImageView : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetImageView(QWidget *parent = nullptr);

  bool LoadImage(const QString& filename, const QString& mask = "",
                 const QList<QPoint>& points = QList<QPoint>(),
                 const QList<RECT_REGION>& regions = QList<RECT_REGION>());

  bool LoadImage(const QString& filename, const QStringList& preprocessed_masks);

  void SetEditedPoints(const QList<QPoint>& pts);

  void paintEvent(QPaintEvent* e);
  void mousePressEvent(QMouseEvent* e);
  void mouseMoveEvent(QMouseEvent* e);
  void mouseReleaseEvent(QMouseEvent* e);
  void wheelEvent(QWheelEvent* e);
  void resizeEvent(QResizeEvent* e);

  bool eventFilter(QObject *watched, QEvent *event);

  QList<QPoint> GetEditedPoints()
  {
    return m_listPoints;
  }

  QList<RECT_REGION> GetEditedRegions()
  {
    return m_listRegions;
  }

  QString GetFilename()
  {
    return m_sFilename;
  }

  QString GetMaskFilename()
  {
    return m_sMaskFilename;
  }

  QImage ReadImageWithExifAwareness(const QString& filename);

  QList<QImage> GetSelectedMasks()
  {
    return m_listSelectedMasks;
  }

  double GetMaskOpacity()
  {
    return m_dMaskOpacity;
  }

  enum EditMode { EM_POINT = 0, EM_REGION, EM_CALIBRATION, EM_SELECT_MASK, EM_EDIT_MASK };

signals:
  void LastRegionEdited(int n);
  void CalibrationReady(const QList<QPoint>& pts);

public slots:
  void SetNumberOfExpectedPoints(int n)
  {
    m_nNumberOfExpectedPoints = n;
  }

  void SetEditMode(int n);

  void ShowMessage(const QString& msg)
  {
    m_sMessage = msg;
    update();
  }

  void HideMessage()
  {
    m_sMessage.clear();
    update();
  }

  void SetOverlay(const QImage& overlay_image);
  void ClearOverlay()
  {
    SetOverlay(QImage());
  }

  void ClearEdits();

  void SetMaskOpacity(double val);

  void SetBrushSize(int n)
  {
    m_nBrushSize = n;
  }

  void SaveMaskIfEdited();

  void UndoLastMaskEdit();

private:
  void PrepareImage();
  void UpdateScaledImage(bool bSmooth = false);
  QPoint ScreenToImage(const QPoint& pt);
  void SetAlphaByMask(QImage& image);
  void FreeHandOnMaskImage(const QPoint& scr_pt1, const QPoint& scr_pt2);
  void UpdatePointOnMaskImage(const QPoint& pt_in);
  void FloodFillMaskImage(const QPoint& scr_pt);
  void UpdateAll();

  QString   m_sFilename;
  QString   m_sMaskFilename;
  QImage    m_image;
  QImage    m_imageMask;
  QImage    m_imageScaled;
  QImage    m_imageOverlay;
  QImage    m_imageOriginal;
  double    m_dScale;
  double    m_dOldScale;
  QPoint    m_ptOffset;
  QPoint    m_ptOldOffset;
  QPoint    m_ptPress;
  QPoint    m_ptPrev;
  bool      m_bPanning;
  bool      m_bZooming;
  bool      m_bDrawing;
  bool      m_bErasing;

  int       m_nNumberOfExpectedPoints;
  int       m_nEditMode;
  QList<QPoint> m_listPoints;
  QList<RECT_REGION> m_listRegions;
  QColor    m_colorPen;

  QString   m_sMessage;
  QStringList m_listPreMasks;
  QImage    m_imageCombinedMaskOverlay;
  double    m_dMaskOpacity;
  QList<QImage> m_listAllMasks;
  QList<QImage> m_listSelectedMasks;
  int      m_nBrushSize;

  bool    m_bMaskEdited;
  QList<QImage>  m_listMaskUndoBuffer;
};

#endif // WIDGETIMAGEVIEW_H
