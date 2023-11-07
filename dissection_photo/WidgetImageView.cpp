#include "WidgetImageView.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "math.h"
#include "exif.h"
#include <QFile>

WidgetImageView::WidgetImageView(QWidget *parent)
  : QWidget(parent), m_dScale(1.0), m_ptOffset(QPoint(0,0)), m_bPanning(false), m_bZooming(false), m_bDrawing(false),
    m_nNumberOfExpectedPoints(2), m_nEditMode(0), m_dMaskOpacity(0.7)
{
  setMouseTracking(true);
  m_colorPen = QColor(50,255,50);

  if (parent)
    parent->installEventFilter(this);
}

bool WidgetImageView::eventFilter(QObject *watched, QEvent *event)
{
  if (watched == parentWidget() && event->type() == QEvent::KeyRelease)
  {
    QKeyEvent* e = static_cast<QKeyEvent*>(event);
    if (e->key() == Qt::Key_Control)
    {
      unsetCursor();
    }
  }

  return QWidget::eventFilter(watched, event);
}

QImage WidgetImageView::ReadImageWithExifAwareness(const QString& filename)
{
  QImage image(filename);
  if (!image.isNull())
  {
    QFile file(filename);
    if (file.open(QIODevice::ReadOnly))
    {
      Exif exif;
      int orientation = 0;
      exif.getExifOrientation(file, &orientation);
      file.close();
      QTransform t;
      if (orientation == 7 || orientation == 8)
        t.rotate(-90);
      else if (orientation == 3 || orientation == 4)
        t.rotate(180);
      else if (orientation == 5 || orientation == 6)
        t.rotate(90);
      if (orientation > 1)
        image = image.transformed(t, Qt::SmoothTransformation);
    }
  }
  return image;
}

bool WidgetImageView::LoadImage(const QString& filename, const QString& mask, const QList<QPoint>& points, const QList<RECT_REGION>& regions)
{
  m_sMessage.clear();
  m_sFilename = filename;
  m_sMaskFilename = mask;
  m_imageOverlay = QImage();
  PrepareImage();
  if (!m_image.isNull())
  {
    m_dScale = 1.0;
    m_ptOffset = QPoint(0,0);
    m_listPoints = points;
    m_listRegions = regions;
    UpdateScaledImage();
    return true;
  }
  else
    return false;
}

bool WidgetImageView::LoadImage(const QString &filename, const QStringList &preprocessed_masks)
{
  m_sMessage.clear();
  m_sFilename = filename;
  m_listPreMasks = preprocessed_masks;
  m_imageCombinedMaskOverlay = QImage();
  PrepareImage();
  if (!m_image.isNull())
  {
    m_dScale = 1.0;
    m_ptOffset = QPoint(0,0);
    UpdateScaledImage();
    return true;
  }
  else
    return false;
}

void WidgetImageView::SetAlphaByMask(QImage& image)
{
  for (int y = 0; y < image.height(); y++)
  {
    QRgb* p_dest = (QRgb*)image.scanLine(y);
    for (int x = 0; x < image.width(); x++)
    {
      if (qRed(p_dest[x]) == 0)
      {
        p_dest[x] = qRgba(0,0,0,0);
      }
    }
  }
}

void WidgetImageView::PrepareImage()
{
  QImage image = ReadImageWithExifAwareness(m_sFilename);
  if (!image.isNull())
  {
    m_image = image;
    if (!m_sMaskFilename.isEmpty())
    {
      QImage mask_image = ReadImageWithExifAwareness(m_sMaskFilename);
      QPainter p(&m_image);
      p.setCompositionMode(QPainter::CompositionMode_Multiply);
      p.drawImage(0, 0, mask_image);
      p.end();
    }
    if (!m_imageOverlay.isNull())
    {
      QPainter p(&m_image);
      p.drawImage(0, 0, m_imageOverlay);
      p.end();
    }
    if (m_imageCombinedMaskOverlay.isNull() && !m_listPreMasks.isEmpty())
    {
      QList<QColor> colors;
      colors << QColor(255,100,100) << QColor(255,255,100) << QColor(100,255,100)
             << QColor(110,245,255) << QColor(75,100,255) << QColor(255,128,0)
             << QColor(100,150,170) << QColor(120,60,220);
      m_imageCombinedMaskOverlay = QImage(m_image.size(), QImage::Format_ARGB32);
      m_imageCombinedMaskOverlay.fill(QColor(0,0,0,0));
      m_listAllMasks.clear();
      m_listSelectedMasks.clear();
      QPainter p;
      p.begin(&m_imageCombinedMaskOverlay);
      for (int i = 0; i < m_listPreMasks.size(); i++)
      {
        QColor c = colors[i%(colors.size())];
        QImage img(m_listPreMasks[i]);
        img = img.convertedTo(QImage::Format_ARGB32);
        SetAlphaByMask(img);
        m_listAllMasks << img;
        QPainter p1(&img);
        p1.setCompositionMode(QPainter::CompositionMode_SourceIn);
        p1.fillRect(img.rect(), c);
        p1.end();
        p.drawImage(0, 0, img);
      }
      p.end();
    }
    if (!m_imageCombinedMaskOverlay.isNull())
    {
      QPainter p(&m_image);
      p.setOpacity(m_dMaskOpacity);
      p.drawImage(0, 0, m_imageCombinedMaskOverlay);
      p.end();
    }
    if (!m_listSelectedMasks.isEmpty())
    {
      QPainter p(&m_image);
      double val = 1;
      p.setOpacity(qMin(1.0, m_dMaskOpacity*2));
      foreach (QImage img, m_listSelectedMasks)
        p.drawImage(0, 0, img);
      p.end();
    }
  }
}

void WidgetImageView::SetEditMode(int n)
{
  m_nEditMode = n;
  if (n == EM_CALIBRATION)
    m_nNumberOfExpectedPoints = 8;
}

void WidgetImageView::SetOverlay(const QImage& overlay_image)
{
  m_imageOverlay = overlay_image;
  PrepareImage();
  UpdateScaledImage();
}

void WidgetImageView::paintEvent(QPaintEvent *e)
{
  QRect rc = rect();
  QPainter p(this);
  p.fillRect(rc, Qt::black);
  QRect target = m_imageScaled.rect();
  target.moveCenter(rect().center()+m_ptOffset);
  p.drawImage(target, m_imageScaled);
  if (m_nEditMode == EM_POINT)
  {
    p.setPen(m_colorPen);
    p.setBrush(m_colorPen);
    QList<QPoint> pts;
    foreach (QPoint pt, m_listPoints)
    {
      pt = pt*m_imageScaled.width()/m_image.width() + target.topLeft();
      p.drawEllipse(pt, 2, 2);
      pts << pt;
    }
    if (pts.size() > 1 && m_nNumberOfExpectedPoints == 2)
    {
      p.setPen(QPen(m_colorPen, 2));
      p.drawLine(pts[0], pts[1]);
    }
  }
  else if (m_nEditMode == EM_REGION)
  {
    QColor pen_color(255,255,255);
    QColor active_color(255,30,30);
    p.setPen(QPen(pen_color, 2));
    p.setBrush(Qt::NoBrush);
    for (int i = 0; i < m_listRegions.size(); i++)
    {
      if (i == m_listRegions.size()-1 && m_bDrawing)
        p.setPen(QPen(active_color, 2));

      QRect rc(m_listRegions[i].first*m_imageScaled.width()/m_image.width() + target.topLeft(),
               m_listRegions[i].second*m_imageScaled.width()/m_image.width() + target.topLeft());
      p.drawRect(rc);
    }
  }
  else if (m_nEditMode == EM_CALIBRATION)
  {
    for (int i = 0; i < m_listPoints.size(); i++)
    {
      QPoint pt = m_listPoints[i];
      pt = pt*m_imageScaled.width()/m_image.width() + target.topLeft();
      if (i%2 == 0)
      {
        p.setBrush(m_colorPen);
        p.setPen(m_colorPen);
        p.drawEllipse(pt, 2, 2);
      }
      else if (i%2 == 1)
      {
        QPoint prev_pt = m_listPoints[i-1];
        prev_pt = prev_pt*m_imageScaled.width()/m_image.width() + target.topLeft();
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(m_colorPen,2));
        double r = sqrt((pt.x()-prev_pt.x())*(pt.x()-prev_pt.x())+(pt.y()-prev_pt.y())*(pt.y()-prev_pt.y()));
        p.drawEllipse(QPointF(prev_pt), r, r);
      }
    }
  }
  if (!m_sMessage.isEmpty())
  {
    QFont f = p.font();
    f.setPixelSize(16);
    p.setFont(f);
    QRect rc = p.fontMetrics().boundingRect(m_sMessage);
    rc.adjust(-20,-9, 20, 9);
    rc.moveCenter(rect().center());
    rc.moveBottom(rect().bottom()-18);
    p.setPen(Qt::NoPen);
    p.setBrush(QColor(0,0,0,200));
    p.drawRoundedRect(rc, rc.height()/2, rc.height()/2);
    p.setPen(Qt::white);
    p.drawText(rc, m_sMessage, Qt::AlignHCenter|Qt::AlignVCenter);
  }
}

void WidgetImageView::mousePressEvent(QMouseEvent *e)
{
  m_ptPress = e->pos();
  m_ptOldOffset = m_ptOffset;
  m_dOldScale = m_dScale;
  if (e->button() == Qt::LeftButton)
  {
    if (e->modifiers() & Qt::ControlModifier)
    {
      m_bDrawing = true;
      if (m_nEditMode != EM_SELECT_MASK)
      {
        QPoint pt = ScreenToImage(m_ptPress);
        RECT_REGION region;
        region.first = pt;
        region.second = pt;
        m_listRegions << region;
        setCursor(Qt::CrossCursor);
      }
    }
    else
      m_bPanning = true;
  }
  else if (e->button() == Qt::RightButton)
    m_bZooming = true;
}

QPoint WidgetImageView::ScreenToImage(const QPoint &pt_in)
{
  QRect target = m_imageScaled.rect();
  target.moveCenter(rect().center()+m_ptOffset);
  QPoint pt = (pt_in - target.topLeft());
  pt = pt*m_image.width()/m_imageScaled.width();
  return pt;
}

void WidgetImageView::mouseReleaseEvent(QMouseEvent *e)
{
  QPoint dpt = e->pos() - m_ptPress;
  if (m_bDrawing && !m_imageScaled.isNull())
  {
    if ((m_nEditMode == EM_POINT || m_nEditMode == EM_CALIBRATION) && qAbs(dpt.x()) < 2 && qAbs(dpt.y()) < 2)
    {
      QPoint pt = ScreenToImage(m_ptPress);
      if (m_listPoints.size() < m_nNumberOfExpectedPoints)
      {
        m_listPoints << pt;
      }
      else
        m_listPoints[m_nNumberOfExpectedPoints-1] = pt;

      if (m_listPoints.size() == m_nNumberOfExpectedPoints && m_nEditMode == EM_CALIBRATION)
        emit CalibrationReady(m_listPoints);

      update();
    }
    else if (m_nEditMode == EM_REGION)
    {
      update();
      QPoint pt = m_listRegions.last().first - m_listRegions.last().second;
      if (qAbs(pt.x()) < 3 || qAbs(pt.y()) < 3)
        m_listRegions.removeLast();
      else
        emit LastRegionEdited(m_listRegions.size()-1);
    }
    else if (m_nEditMode == EM_SELECT_MASK && qAbs(dpt.x()) < 3 && qAbs(dpt.y() < 3))
    {
      QPoint pt = ScreenToImage(m_ptPress);
      for (int i = m_listAllMasks.size()-1; i >= 0; i--)
      {
        QImage img = m_listSelectedMasks[i];
        if (qRed(img.pixel(pt.x(), pt.y())) > 0 && !m_listSelectedMasks.contains(img))
        {
          m_listSelectedMasks << img;
          PrepareImage();
          UpdateScaledImage();
          update();
          break;
        }
      }
    }
  }
  else if (m_bZooming)
  {
    UpdateScaledImage(true);
  }
  m_bPanning = false;
  m_bZooming = false;
  m_bDrawing = false;
  unsetCursor();
}

void WidgetImageView::mouseMoveEvent(QMouseEvent *e)
{
  QPoint dpt = e->pos() - m_ptPress;
  if (m_bPanning)
  {
    m_ptOffset = m_ptOldOffset + dpt;
    update();
  }
  else if (m_bZooming)
  {
    int dn = e->y()-m_ptPress.y();
    if (dn > 0)
      m_dScale = m_dOldScale/pow(1.002, dn);
    else
      m_dScale = m_dOldScale*pow(1.002, -dn);
    if (m_dScale < 0.25)
      m_dScale = 0.25;
    else if (m_dScale > 10)
      m_dScale = 8;
    UpdateScaledImage();
  }
  else if (m_bDrawing && m_nEditMode == EM_REGION && !m_listRegions.isEmpty())
  {
    QPoint pt_0 = ScreenToImage(m_ptPress), pt_1 = ScreenToImage(e->pos());
    if (pt_0.x() > pt_1.x())
    {
      int temp = pt_0.x();
      pt_0.setX(pt_1.x());
      pt_1.setX(temp);
    }
    if (pt_0.y() > pt_1.y())
    {
      int temp = pt_0.y();
      pt_0.setY(pt_1.y());
      pt_1.setY(temp);
    }
    RECT_REGION region = m_listRegions.last();
    region.first = pt_0;
    region.second = pt_1;
    m_listRegions[m_listRegions.size()-1] = region;
    update();
  }
  if (e->modifiers() & Qt::ControlModifier)
  {
    setCursor(Qt::CrossCursor);
  }

  QWidget::mouseMoveEvent(e);
}

void WidgetImageView::wheelEvent(QWheelEvent* e)
{
  //  qDebug() << e;
}

void WidgetImageView::resizeEvent(QResizeEvent *e)
{
  UpdateScaledImage();
}

void WidgetImageView::UpdateScaledImage(bool bSmooth)
{
  if (!m_image.isNull() && height() > 0)
  {
    if (1.0*m_image.width()/m_image.height() > 1.0*width()/height())
      m_imageScaled = m_image.scaledToWidth(width()*m_dScale, bSmooth?Qt::SmoothTransformation:Qt::FastTransformation);
    else
      m_imageScaled = m_image.scaledToHeight(height()*m_dScale, bSmooth?Qt::SmoothTransformation:Qt::FastTransformation);
  }
  update();
}

void WidgetImageView::ClearEdits()
{
  m_listPoints.clear();
  m_listRegions.clear();
  m_imageOverlay = QImage();
  m_listSelectedMasks.clear();
  PrepareImage();
  UpdateScaledImage();
  HideMessage();
  update();
}

void WidgetImageView::SetMaskOpacity(double val)
{
  m_dMaskOpacity = val;
  PrepareImage();
  UpdateScaledImage();
  update();
}
