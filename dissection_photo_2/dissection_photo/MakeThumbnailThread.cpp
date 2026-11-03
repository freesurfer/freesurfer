#include "MakeThumbnailThread.h"
#include <QPainter>
#include <QFontMetrics>
#include "MyUtils.h"
#include "MaskProcessor.h"

MakeThumbnailThread::MakeThumbnailThread(QObject *parent) : QThread(parent)
{
  m_ratio = 1;
}

void MakeThumbnailThread::run()
{
  m_bAbort = false;
  if (m_nUpdateIndex >= 0)
  {
    QPixmap pix = !m_strUpdateLabels.isEmpty() ? BuildComponentsPixmap(m_nUpdateIndex)
                : m_strUpdateFile.isEmpty() ? BuildPixmap(m_nUpdateIndex, m_images.first())
                                            : BuildPixmapFromFiles(m_nUpdateIndex, m_strUpdateFile, m_strUpdateMask);
    emit PixmapReady(pix, m_nUpdateIndex, m_nGeneration);
  }
  else
  {
    int ncnt = m_files.size();
    for (int i = 0; i < ncnt; i++)
    {
      QPixmap pix = BuildPixmapFromFiles(i, m_files[i], i < m_maskFiles.size() ? m_maskFiles[i] : QString());
      if (m_bAbort)
        break;
      emit PixmapReady(pix, -1, m_nGeneration);
    }
  }
}

QPixmap MakeThumbnailThread::BuildPixmapFromFiles(int nIndex, const QString& file, const QString& mask_file)
{
  // The mask (NpyToImage()) has no EXIF tag of its own: give it the photo's
  // orientation, or it would stay unrotated while the photo rotates.
  const int nOrientation = MyUtils::GetExifOrientation(file);
  const QImage img = MyUtils::ApplyExifOrientation(QImage(file), nOrientation);
  QImage mask_img;
  if (!mask_file.isEmpty())
  {
    mask_img = MyUtils::ApplyExifOrientation(QImage(mask_file), nOrientation).convertToFormat(QImage::Format_ARGB32);
    for (int y = 0; y < mask_img.height(); y++)
    {
      QRgb* p = (QRgb*)mask_img.scanLine(y);
      for (int x = 0; x < mask_img.width(); x++)
      {
        if (qRed(p[x]) == 0)
        {
          p[x] = qRgba(0,0,0,0);
        }
      }
    }
  }
  return BuildPixmap(nIndex, img, mask_img);
}

QPixmap MakeThumbnailThread::BuildPixmap(int nIndex, const QImage &img_in, const QImage& mask_in)
{
  // A null or empty source would divide by zero and paint on a null pixmap:
  // return a blank placeholder instead.
  if (img_in.isNull() || img_in.width() <= 0 || img_in.height() <= 0)
  {
    QPixmap pix(m_size);
    pix.fill(Qt::darkGray);
    pix.setDevicePixelRatio(m_ratio);
    return pix;
  }

  QImage img = img_in;
  if (!mask_in.isNull())
  {
    QPainter p(&img);
    QImage mask = mask_in;
    // Defensive: masks normally match the photo size already. drawImage()
    // blits at the mask's own size, so a larger stale mask would cover the
    // whole thumbnail. IgnoreAspectRatio: the mask must land on the photo's
    // pixel grid.
    if (mask.size() != img.size())
      mask = mask.scaled(img.size(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    QPainter p2(&mask);
    p2.setCompositionMode(QPainter::CompositionMode_SourceIn);
    p2.fillRect(mask.rect(), QColor(0,255,0));
    p2.end();

    // 0.3, lighter than the edit view's 0.7: masks usually cover most of the
    // frame, and at 0.7 small thumbnails become solid green.
    p.setOpacity(0.3);
    p.drawImage(0, 0, mask);
    p.end();
  }
  QSize sz = m_size;
  QRect rc = img.rect();
  if (img.width()*1.0/img.height() > sz.width()*1.0/sz.height())
  {
    rc.setWidth(sz.width()*img.height()/sz.height());
  }
  else
  {
    rc.setHeight(sz.height()*img.width()/sz.width());
  }
  rc.moveCenter(img.rect().center());
  img = img.copy(rc);
  QPixmap pix = QPixmap::fromImage(img.scaledToWidth(sz.width()*m_ratio, Qt::SmoothTransformation));
  QPainter p(&pix);

  // 1-based image number, drawn as a pill badge in the corner.
  // Everything is scaled by m_ratio: the DPR is set on the pixmap only after
  // painting, so coordinates here are device pixels.
  p.setRenderHint(QPainter::Antialiasing, true);

  const QString strIndex = QString::number(nIndex+1);
  QFont f = p.font();
  // setPixelSize: point sizes depend on the paint device's DPI, which isn't
  // reliable for an offscreen pixmap on a worker thread.
  f.setPixelSize(12*m_ratio);
  f.setBold(true);
  p.setFont(f);
  // Metrics from the painter, so they match what is drawn on this pixmap.
  const QFontMetricsF fm(p.fontMetrics());

  const qreal dPen   = 1.0*m_ratio;   // hairline outline
  const qreal dInset = 3.0*m_ratio;   // gap from the thumbnail's top-left edge
  const qreal dPadX  = 4.5*m_ratio;   // breathing room either side of the digits
  // Fixed pill height, grown only if the font is taller (drawText() clips).
  // Kept small: carousel cells are only 68px.
  const qreal dH     = qMax(17.0*m_ratio, fm.height() + 3.0*m_ratio);
  // Width follows the digits; the qMax keeps a single digit a circle.
  const qreal dW     = qMax(dH, fm.horizontalAdvance(strIndex) + 2.0*dPadX);
  const QRectF rcBadge(dInset, dInset, dW, dH);

  // A soft drop shadow one pixel below, so the badge separates from the photo.
  p.setPen(Qt::NoPen);
  p.setBrush(QColor(0, 0, 0, 70));
  p.drawRoundedRect(rcBadge.translated(0, dPen), dH/2.0, dH/2.0);

  // Translucent dark body with white digits reads on light and dark photos; the
  // light hairline outline keeps it visible on black.
  const QRectF rcBody = rcBadge.adjusted(dPen/2.0, dPen/2.0, -dPen/2.0, -dPen/2.0);
  p.setBrush(QColor(0, 0, 0, 170));
  p.setPen(QPen(QColor(255, 255, 255, 115), dPen));
  p.drawRoundedRect(rcBody, rcBody.height()/2.0, rcBody.height()/2.0);

  p.setBrush(Qt::NoBrush);
  p.setPen(Qt::white);
  p.drawText(rcBadge, Qt::AlignCenter, strIndex);
  p.end();
  pix.setDevicePixelRatio(m_ratio);
  return pix;
}

void MakeThumbnailThread::Abort()
{
  m_mutex.lock();
  m_bAbort = true;
  m_mutex.unlock();
}

void MakeThumbnailThread::UpdateImage(int n, const QImage &img)
{
  m_images.clear();
  m_images << img;
  m_strUpdateFile.clear();
  m_strUpdateMask.clear();
  m_strUpdateLabels.clear();
  m_nUpdateIndex = n;
  start();
}

void MakeThumbnailThread::UpdateFiles(int n, const QString& file, const QString& mask_file)
{
  m_images.clear();
  m_strUpdateFile = file;
  m_strUpdateMask = mask_file;
  m_strUpdateLabels.clear();
  m_nUpdateIndex = n;
  start();
}

void MakeThumbnailThread::UpdateComponents(int n, const QString& file, const QString& mask_file, const QString& labels_file,
                                           const QList<QColor>& colors, double dOpacity)
{
  m_images.clear();
  m_strUpdateFile = file;
  m_strUpdateMask = mask_file;
  m_strUpdateLabels = labels_file;
  m_listColors = colors;
  m_dOpacity = dOpacity;
  m_nUpdateIndex = n;
  start();
}

QPixmap MakeThumbnailThread::BuildComponentsPixmap(int nIndex)
{
  // As WidgetImageView::PrepareImage() composes the Connected components view:
  // the photo multiplied by its mask, then the components at the overlay opacity.
  QImage img = MyUtils::ReadImageWithExifAwareness(m_strUpdateFile).convertToFormat(QImage::Format_ARGB32);
  if (img.isNull())
    return BuildPixmap(nIndex, img);
  const QImage mask = MyUtils::ApplyExifOrientation(QImage(m_strUpdateMask), MyUtils::GetExifOrientation(m_strUpdateFile))
                          .convertToFormat(QImage::Format_ARGB32);
  QPainter p(&img);
  if (!mask.isNull())
  {
    p.setCompositionMode(QPainter::CompositionMode_Multiply);
    p.drawImage(0, 0, mask);
    p.setCompositionMode(QPainter::CompositionMode_SourceOver);
  }
  const QImage labels = MaskProcessor::SavedMaskImage(m_strUpdateLabels, m_listColors);
  if (!labels.isNull())
  {
    p.setOpacity(m_dOpacity);
    p.drawImage(0, 0, labels);
  }
  p.end();
  return BuildPixmap(nIndex, img);
}
