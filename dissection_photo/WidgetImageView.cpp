#include "WidgetImageView.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "math.h"
#include "exif.h"
#include <QFile>
#include "MyUtils.h"
#include <QGuiApplication>

WidgetImageView::WidgetImageView(QWidget *parent)
  : QWidget(parent), m_dScale(1.0), m_ptOffset(QPoint(0,0)), m_bPanning(false), m_bZooming(false), m_bDrawing(false),
    m_nNumberOfExpectedPoints(2), m_nEditMode(0), m_dMaskOpacity(0.7), m_nBrushSize(1), m_bMaskEdited(false)
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
  m_listMaskUndoBuffer.clear();
  m_sFilename = filename;
  m_sMaskFilename = mask;
  m_imageOriginal = m_imageMask = QImage();
  if (!mask.isEmpty())
    m_imageMask = ReadImageWithExifAwareness(m_sMaskFilename).convertToFormat(QImage::Format_ARGB32);
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
  m_listMaskUndoBuffer.clear();
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
  if (m_imageOriginal.isNull())
    m_imageOriginal = ReadImageWithExifAwareness(m_sFilename);
  if (!m_imageOriginal.isNull())
  {
    m_image = m_imageOriginal;
    QImage image = m_imageMask;
    if (!m_sMaskFilename.isEmpty())
    {
      QPainter p(&m_image);
      if (m_nEditMode == EM_EDIT_MASK)
      {
        p.setOpacity(m_dMaskOpacity);
        QPainter p2(&image);
        p2.setCompositionMode(QPainter::CompositionMode_Darken);
        p2.fillRect(image.rect(), QColor(0,255,0));
      }
      else
      {
        p.setCompositionMode(QPainter::CompositionMode_Multiply);
      }
      p.drawImage(0, 0, image);
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
        img = img.convertToFormat(QImage::Format_ARGB32);
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
  else if (m_nEditMode == EM_EDIT_MASK)
  {
    if (QGuiApplication::queryKeyboardModifiers() & Qt::ControlModifier)
    {
      QPointF pt = this->mapFromGlobal(QCursor::pos());
      p.setBrush(QColor(255,0,0, 100));
      p.setPen(Qt::NoPen);
      double r = m_nBrushSize/2.0*m_imageScaled.width()/m_image.width();
      p.drawEllipse(pt, r, r);
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
  m_ptPrev = e->pos();
  m_ptOldOffset = m_ptOffset;
  m_dOldScale = m_dScale;
  if (e->button() == Qt::LeftButton)
  {
    if (e->modifiers() & Qt::ControlModifier)
    {
      m_bDrawing = true;
      if (m_nEditMode == EM_EDIT_MASK)
      {
        m_bErasing = false;
        m_listMaskUndoBuffer << m_imageMask;
        FreeHandOnMaskImage(m_ptPrev, m_ptPrev);
      }
      else if (m_nEditMode != EM_SELECT_MASK)
      {
        QPoint pt = ScreenToImage(m_ptPress);
        RECT_REGION region;
        region.first = pt;
        region.second = pt;
        m_listRegions << region;
        setCursor(Qt::CrossCursor);
      }
    }
    else if (e->modifiers() & Qt::ShiftModifier)
    {
      m_bDrawing = true;
      if (m_nEditMode == EM_EDIT_MASK)
      {
        m_listMaskUndoBuffer << m_imageMask;
        FloodFillMaskImage(m_ptPress);
      }
    }
    else
      m_bPanning = true;
  }
  else if (e->button() == Qt::RightButton)
  {
    if (e->modifiers() & Qt::ControlModifier)
    {
      if (m_nEditMode == EM_EDIT_MASK)
      {
        m_bDrawing = true;
        m_bErasing = true;
        m_listMaskUndoBuffer << m_imageMask;
        FreeHandOnMaskImage(m_ptPrev, m_ptPrev);
      }
    }
    else
      m_bZooming = true;
  }
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
          UpdateAll();
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
  m_bErasing = false;
  unsetCursor();
  update();
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
  else if (m_bDrawing && m_nEditMode == EM_EDIT_MASK)
  {
    FreeHandOnMaskImage(m_ptPrev, e->pos());
    m_ptPrev = e->pos();
  }
  if (e->modifiers() & Qt::ControlModifier)
  {
    if (m_nEditMode == EM_EDIT_MASK)
      update();
    else
      setCursor(Qt::CrossCursor);
  }

  QWidget::mouseMoveEvent(e);
}

void WidgetImageView::keyPressEvent(QKeyEvent* e)
{
  if (m_bDrawing && m_nEditMode == EM_EDIT_MASK && e->modifiers() & Qt::ControlModifier)
  {
    update();
  }
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
  HideMessage();
  UpdateAll();
}

void WidgetImageView::SetMaskOpacity(double val)
{
  m_dMaskOpacity = val;
  UpdateAll();
}

void WidgetImageView::SetEditedPoints(const QList<QPoint>& pts)
{
  m_listPoints = pts;
  m_imageOverlay = QImage();
  m_listSelectedMasks.clear();
  UpdateAll();
  HideMessage();
}

void WidgetImageView::FreeHandOnMaskImage(const QPoint& scr_pt1, const QPoint& scr_pt2)
{
  QPoint pt1 = ScreenToImage(scr_pt1), pt2 = ScreenToImage(scr_pt2);
  if (pt1 == pt2)
    UpdatePointOnMaskImage(pt1);
  else
  {
    int x0 = pt1.x(), x1 = pt2.x(), y0 = pt1.y(), y1 = pt2.y();
    int dx = x1 - x0;
    int dy = y1 - y0;
    double t = 0.5;
//    list = SetVoxelByIndex( n1, nPlane, bAdd, ignore_brush_size );
    if ( abs( dx ) > abs( dy ) )
    {
      double m = (double) dy / (double) dx;
      t += y0;
      dx = ( dx < 0 ? -1 : 1 );
      m *= dx;
      while ( x0 != x1 )
      {
        x0 += dx;
        t += m;
        UpdatePointOnMaskImage(QPoint(x0, t));
      }
    }
    else
    {
      double m = (double) dx / (double) dy;
      t += x0;
      dy = ( dy < 0 ? -1 : 1 );
      m *= dy;
      while ( y0 != y1 )
      {
        y0 += dy;
        t += m;
        UpdatePointOnMaskImage(QPoint(t, y0));
      }
    }
  }
  UpdateAll();
}

void WidgetImageView::UpdatePointOnMaskImage(const QPoint& pt_in)
{
  int nr = m_nBrushSize/2+1;
  int w = m_imageMask.width(), h = m_imageMask.height();
  QRgb fill_val = QColor(m_bErasing?Qt::black:Qt::white).rgb();
  for (int i = -nr+1; i < nr; i++)
  {
    for (int j = -nr+1; j < nr; j++)
    {
      QPoint pt(pt_in.x()+i, pt_in.y()+j);
      if (pt.x() >= 0 && pt.x() < w && pt.y() >= 0 && pt.y() < h &&
          sqrt(i*i+j*j) <= m_nBrushSize/2.0)
      {
        QRgb* p = (QRgb*)m_imageMask.scanLine(pt.y());
        p[pt.x()] = fill_val;
      }
    }
  }
  m_bMaskEdited = true;
}

void WidgetImageView::SaveMaskIfEdited()
{
  if (m_bMaskEdited)
  {
    m_imageMask.save(m_sMaskFilename);
    m_bMaskEdited = false;
  }
}

void WidgetImageView::FloodFillMaskImage(const QPoint& pt_in)
{
  QPoint pt = ScreenToImage(pt_in);
  int w = m_imageMask.width(), h = m_imageMask.height();
  if (pt.x() < 0 || pt.x() >= w || pt.y() < 0 || pt.y() >= h)
    return;

  char* buf = new char[w*h];
  memset(buf, 0, w*h);
  for (int j = 0; j < h; j++)
  {
    QRgb* p = (QRgb*)m_imageMask.scanLine(j);
    for (int i = 0; i < w; i++)
    {
      if (qRed(p[i]) > 0)
        buf[j*w+i] = 1;
    }
  }
  MyUtils::FloodFill(buf, pt.x(), pt.y(), w, h, 1, 1);
  QRgb fill_val = QColor(Qt::white).rgb();
  for (int j = 0; j < h; j++)
  {
    QRgb* p = (QRgb*)m_imageMask.scanLine(j);
    for (int i = 0; i < w; i++)
    {
      if (buf[j*w+i] > 0)
        p[i] = fill_val;
    }
  }
  delete[] buf;
  m_bMaskEdited = true;
  UpdateAll();
}

void WidgetImageView::UpdateAll()
{
  PrepareImage();
  UpdateScaledImage();
  update();
}

void WidgetImageView::UndoLastMaskEdit()
{
  if (!m_listMaskUndoBuffer.isEmpty())
  {
    m_imageMask = m_listMaskUndoBuffer.last();
    m_listMaskUndoBuffer.removeLast();
    UpdateAll();
  }
}
