#include "QThumbnailWidget.h"
#include <QListWidgetItem>
#include <QDebug>
#include <QPainter>
#include "MakeThumbnailThread.h"

QThumbnailWidget::QThumbnailWidget(QWidget* parent) : QListWidget(parent)
{
  setViewMode(QListWidget::IconMode);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  m_thread = new MakeThumbnailThread;
  connect(m_thread, SIGNAL(PixmapReady(QPixmap,int)), SLOT(OnPixmapReady(QPixmap,int)));

  connect(this, SIGNAL(itemSelectionChanged()), SLOT(OnItemSelectionChanged()), Qt::QueuedConnection);
}

QThumbnailWidget::~QThumbnailWidget()
{
  m_thread->wait();
  m_thread->deleteLater();
}

void QThumbnailWidget::LoadImages(const QList<QImage> &images)
{
  Clear();
  m_thread->SetPixmapSize(iconSize(), devicePixelRatio());
  m_thread->LoadImages(images);
}

void QThumbnailWidget::LoadImageFiles(const QStringList &files, const QStringList& mask_files)
{
  Clear();
  m_thread->SetPixmapSize(iconSize(), devicePixelRatio());
  m_thread->LoadImageFiles(files, mask_files);
}

void QThumbnailWidget::Clear()
{
  if (m_thread->isRunning())
  {
    m_thread->blockSignals(true);
    m_thread->Abort();
    m_thread->wait();
    m_thread->blockSignals(false);
  }
  clear();
  m_listPixmap.clear();
}

void QThumbnailWidget::OnPixmapReady(const QPixmap &pix, int n)
{
  QIcon icn(pix);
  if (n >= 0)
  {
    m_listPixmap[n] = pix;
    QListWidgetItem* it = item(n);
    it->setIcon(icn);
  }
  else
  {
    m_listPixmap << pix;
    QListWidgetItem* item = new QListWidgetItem(icn, "");
    addItem(item);
    if (m_listPixmap.size() == 1)
      SetCurrentIndex(0);
  }
}

void QThumbnailWidget::OnItemSelectionChanged()
{
  emit CurrentIndexChanged(currentRow());
}

void QThumbnailWidget::SetCurrentIndex(int n)
{
  setCurrentRow(n);
}

void QThumbnailWidget::SetItemChecked(int n, bool bChecked)
{
  if (n < 0 || n >= m_listPixmap.size())
    return;

  QPixmap pix = m_listPixmap[n];
  if (bChecked)
  {
    pix.setDevicePixelRatio(1);
    QPainter p(&pix);
    QPixmap pch(":/resource/checked@2x.png");
    QRect rc = pch.rect();
    rc.moveTopRight(pix.rect().topRight()+QPoint(-2,2));
    p.drawPixmap(rc, pch);
    p.end();
    QListWidgetItem* it = item(n);
    it->setIcon(QIcon(pix));
  }
}

void QThumbnailWidget::SetItemImage(int n, const QImage& img)
{
  m_thread->UpdateImage(n, img);
}
