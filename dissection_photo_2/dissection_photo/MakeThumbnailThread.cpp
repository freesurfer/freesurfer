#include "MakeThumbnailThread.h"
#include <QPainter>
#include <QDebug>
#include "MyUtils.h"

MakeThumbnailThread::MakeThumbnailThread(QObject *parent) : QThread(parent)
{
  m_ratio = 1;
}

void MakeThumbnailThread::run()
{
  m_bAbort = false;
  if (m_nUpdateIndex >= 0)
  {
    QPixmap pix = BuildPixmap(m_nUpdateIndex, m_images.first());
    emit PixmapReady(pix, m_nUpdateIndex);
  }
  else
  {
    int ncnt = qMax(m_images.size(), m_files.size());
    for (int i = 0; i < ncnt; i++)
    {
      QImage img, mask_img;
      if (!m_files.isEmpty())
        img = MyUtils::ReadImageWithExifAwareness(m_files[i]);
      else
        img = m_images[i];
      if (i < m_maskFiles.size())
      {
        mask_img = MyUtils::ReadImageWithExifAwareness(m_maskFiles[i]).convertToFormat(QImage::Format_ARGB32);
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
      QPixmap pix = BuildPixmap(i, img, mask_img);
      if (m_bAbort)
        break;
      else
        emit PixmapReady(pix);
    }
  }
}

QPixmap MakeThumbnailThread::BuildPixmap(int nIndex, const QImage &img_in, const QImage& mask_in)
{
  QImage img = img_in;
  if (!mask_in.isNull())
  {
    QPainter p(&img);
    QImage mask = mask_in;
    QPainter p2(&mask);
    p2.setCompositionMode(QPainter::CompositionMode_SourceIn);
    p2.fillRect(mask.rect(), QColor(0,255,0));
    p2.end();

    p.setOpacity(0.7);
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
  QFont f = p.font();
  f.setPointSize(12*m_ratio);
  p.setFont(f);
  QRect frc = QFontMetrics(f).boundingRect(QString::number(nIndex+1));
  frc.adjust(-3*m_ratio, -m_ratio, 3*m_ratio, m_ratio);
  if (frc.width() < frc.height())
    frc.setWidth(frc.height());
  frc.moveTopLeft(QPoint(0,0));
  p.fillRect(frc, Qt::white);
  p.setPen(Qt::black);
  p.drawText(frc, Qt::AlignCenter, QString::number(nIndex+1));
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
  m_files.clear();
  m_maskFiles.clear();
  m_nUpdateIndex = n;
  start();
}
