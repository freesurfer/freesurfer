#include "RenderWidget.h"
#include <QDebug>
#include <QPainter>
#include <QFileInfo>
#include <QMouseEvent>

RenderWidget::RenderWidget(QWidget *parent) :
    QWidget(parent),
    m_bAutoResize(true),
    m_nCurrentImageIndex(-1),
    m_nInterval(40)
{
  connect(&m_timerRender, SIGNAL(timeout()), this, SLOT(OnTimer()));
}

int RenderWidget::LoadImages(const QStringList& filenames)
{
  foreach (QString fn, filenames)
  {
    QImage image(fn);
    if (!image.isNull())
    {
      QFileInfo fi(fn);
      image.setText("FileName", fi.fileName());
      image.setText("FullPath", fi.absoluteFilePath());
      m_images << image;
      m_resizedImages << QImage();
    }
    else
      qWarning("%s can not be loaded as image file", qPrintable(fn));
  }
  if (!m_images.isEmpty())
  {
    SetCurrentImageIndex(0);
    qWarning("Loaded %d image files", (int)m_images.size());
  }

  return m_images.size();
}

void RenderWidget::SetCurrentImageIndex(int n)
{
  m_nCurrentImageIndex = n;
  repaint();
  emit CurrentImageChanged(m_images[n], n);
}

void RenderWidget::paintEvent(QPaintEvent *)
{
  QRect rc = rect();
  QPainter p(this);
  p.fillRect(rc, Qt::black);
  if (m_bAutoResize)
  {
    QImage image = m_resizedImages[m_nCurrentImageIndex];
    if (image.isNull())
    {
      image = m_images[m_nCurrentImageIndex].scaled(rc.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation );
      m_resizedImages[m_nCurrentImageIndex] = image;
    }

    p.drawImage((rc.width()-image.width())/2, (rc.height()-image.height())/2, image);
  }
  else
    p.drawImage(0, 0, m_images[m_nCurrentImageIndex]);
}

void RenderWidget::resizeEvent(QResizeEvent *)
{
  m_resizedImages.clear();
  for (int i = 0; i < m_images.size(); i++)
    m_resizedImages << QImage();
}

void RenderWidget::mousePressEvent(QMouseEvent *e)
{
  if (e->button() == Qt::LeftButton)
  {
    m_bPressed = true;
    m_nY = e->pos().y();
  }
}

void RenderWidget::mouseMoveEvent(QMouseEvent *e)
{
  if (m_bPressed)
  {
    int d = e->pos().y()-m_nY;
    int nStepSize = 1;
    if (d >= nStepSize)
    {
      do {
        OnForward();
        d -= nStepSize;
      }
      while (d > 0);
      m_nY = e->pos().y();
    }
    else if (d <= -nStepSize)
    {
      do {
        OnBack();
        d += nStepSize;
      }
      while (d < 0);
      m_nY = e->pos().y();
    }
  }
}

void RenderWidget::mouseReleaseEvent(QMouseEvent *e)
{
  Q_UNUSED(e);

  m_bPressed = false;
}

void RenderWidget::OnLoop()
{
  m_bSwing = false;
  m_timerRender.start(m_nInterval);
  SetCurrentImageIndex(0);
}

void RenderWidget::OnSwing()
{
  m_bSwing = true;
  SetCurrentImageIndex(0);
  m_timerRender.start(m_nInterval);
}

void RenderWidget::OnBack()
{
  int n = m_nCurrentImageIndex-1;
  if (n < 0)
    n = m_images.size()-1;
  SetCurrentImageIndex(n);
}

void RenderWidget::OnForward()
{
  int n = m_nCurrentImageIndex+1;
  if (n > m_images.size()-1)
    n = 0;
  SetCurrentImageIndex(n);
}

void RenderWidget::OnStop()
{
  m_timerRender.stop();
}

void RenderWidget::OnTimer()
{
  static bool bReverse = false;
  if (m_bSwing)
  {
    int n = 0;
    if (!bReverse)
    {
      n = m_nCurrentImageIndex + 1;
      if (n > m_images.size()-1)
      {
        n = m_images.size()-2;
        bReverse = true;
      }
    }
    else
    {
      n = m_nCurrentImageIndex - 1;
      if (n < 0)
      {
        n = 1;
        bReverse = false;
      }
    }
    SetCurrentImageIndex(n);
  }
  else
  {
    int n = m_nCurrentImageIndex + 1;
    if (n > m_images.size()-1)
      n = 0;
    SetCurrentImageIndex(n);
  }
}

void RenderWidget::SetAutoResize(bool bAuto)
{
  m_bAutoResize = bAuto;
  update();
}

void RenderWidget::SetSpeed(int n)
{
  m_nInterval = 100-n;
  if (m_nInterval < 50)
    m_nInterval++;
  m_timerRender.setInterval(m_nInterval);
}

