#include "QThumbnailWidget.h"
#include <QListWidgetItem>
#include <QPainter>
#include "MakeThumbnailThread.h"

namespace
{
// See ItemPadding(). Applied to both axes so the cell is square.
const int kThumbItemPadding = 12;
}

QThumbnailWidget::QThumbnailWidget(QWidget* parent) : QListWidget(parent)
{
  setViewMode(QListWidget::IconMode);
  // Disable drag and drop. setViewMode(IconMode) silently sets movement to
  // Free, which enables dragging and dropping; an internal drop then moves the
  // item to an arbitrary free spot, stacking thumbnails (looks like duplicates
  // when dragging fast). No view here wants reorderable thumbnails.
  // All three switches are needed:
  //  * setMovement(Static): also stops a later setViewMode(IconMode) (in
  //    ImageGridView or setupUi()) from switching back to Free.
  //  * setDragDropMode(NoDragDrop): clears dragEnabled and the view's
  //    acceptDrops.
  //  * viewport()->setAcceptDrops(false): the viewport is what Free enabled.
  // The folder drop targets (widgetDrop*) are siblings, not ancestors, so they
  // are unaffected.
  setMovement(QListWidget::Static);
  setDragDropMode(QAbstractItemView::NoDragDrop);
  setDragEnabled(false);
  setAcceptDrops(false);
  viewport()->setAcceptDrops(false);
  setDropIndicatorShown(false);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  m_thread = new MakeThumbnailThread;
  connect(m_thread, SIGNAL(PixmapReady(QPixmap,int,int)), SLOT(OnPixmapReady(QPixmap,int,int)));
  connect(m_thread, SIGNAL(finished()), SLOT(OnThreadFinished()), Qt::QueuedConnection);

  connect(this, SIGNAL(itemSelectionChanged()), SLOT(OnItemSelectionChanged()), Qt::QueuedConnection);
}

int QThumbnailWidget::ItemPadding()
{
  return kThumbItemPadding;
}

QThumbnailWidget::~QThumbnailWidget()
{
  m_thread->wait();
  m_thread->deleteLater();
}

void QThumbnailWidget::LoadImageFiles(const QStringList &files, const QStringList& mask_files, const QList<int>& checked_indices)
{
  Clear();
  m_listPendingChecked = checked_indices;
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
  m_listPendingChecked.clear();
  m_mapPendingUpdates.clear();
  m_setChecked.clear();
  // blockSignals() above stops new results, not ones already queued to this thread.
  m_thread->SetGeneration(++m_nGeneration);
}

void QThumbnailWidget::OnPixmapReady(const QPixmap &pix, int n, int nGeneration)
{
  // From a list that has been replaced since.
  if (nGeneration != m_nGeneration)
    return;
  QIcon icn(pix);
  if (n >= 0)
  {
    // An update for an item that is gone (the list was reloaded).
    if (n >= m_listPixmap.size())
      return;
    m_listPixmap[n] = pix;
    item(n)->setIcon(icn);
    if (m_setChecked.contains(n))
      SetItemChecked(n, true);
  }
  else
  {
    m_listPixmap << pix;
    QListWidgetItem* item = new QListWidgetItem(icn, "");
    // Explicit size hint: otherwise IconMode reserves a text line under the
    // empty-text icon (80x99 for a 68x68 icon). Must be set before addItem();
    // later it doesn't change the laid-out rect.
    item->setSizeHint(QSize(iconSize().width() + kThumbItemPadding,
                            iconSize().height() + kThumbItemPadding));
    addItem(item);
    if (m_listPixmap.size() == 1)
      SetCurrentIndex(0);
    // Thumbnails arrive in order, so this one's index is m_listPixmap.size()-1;
    // apply its checked state now.
    if (m_listPendingChecked.contains(m_listPixmap.size()-1) || m_setChecked.contains(m_listPixmap.size()-1))
      SetItemChecked(m_listPixmap.size()-1, true);
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
  // Remembered even before the item has loaded; applied when it arrives.
  if (bChecked)
    m_setChecked.insert(n);
  else
  {
    m_setChecked.remove(n);
    // Also cancels a tick asked for when the list was loaded but not shown yet.
    m_listPendingChecked.removeAll(n);
  }
  if (n < 0 || n >= m_listPixmap.size())
    return;
  // m_listPixmap holds the unmarked thumbnail: unchecking shows it again.
  if (!bChecked)
  {
    item(n)->setIcon(QIcon(m_listPixmap[n]));
    return;
  }

  QPixmap pix = m_listPixmap[n];
  pix.setDevicePixelRatio(1);
  QPainter p(&pix);
  // Light green outline marks checked (processed) images, alongside the
  // checkmark badge.
  p.setRenderHint(QPainter::Antialiasing, true);
  QPen pen(QColor(0x8b, 0xd4, 0x8b));
  int nPenWidth = 3;
  pen.setWidth(nPenWidth);
  p.setPen(pen);
  p.setBrush(Qt::NoBrush);
  p.drawRect(pix.rect().adjusted(nPenWidth/2, nPenWidth/2, -nPenWidth/2-1, -nPenWidth/2-1));

  QPixmap pch(":/resource/checked@2x.png");
  if (devicePixelRatio() == 1)
    pch = pch.scaledToWidth(pch.width()/2, Qt::SmoothTransformation);
  QRect rc = pch.rect();
  int nOffset = devicePixelRatio();
  rc.moveTopRight(pix.rect().topRight()+QPoint(-nOffset,nOffset));
  p.drawPixmap(rc, pch);
  p.end();
  item(n)->setIcon(QIcon(pix));
}

void QThumbnailWidget::SetItemImage(int n, const QImage& img)
{
  Update u;
  u.image = img;
  m_mapPendingUpdates[n] = u;
  StartNextUpdate();
}

void QThumbnailWidget::SetItemFiles(int n, const QString& file, const QString& mask_file)
{
  Update u;
  u.file = file;
  u.mask_file = mask_file;
  m_mapPendingUpdates[n] = u;
  StartNextUpdate();
}

void QThumbnailWidget::SetItemComponents(int n, const QString& file, const QString& mask_file, const QString& labels_file,
                                         const QList<QColor>& colors, double dOpacity)
{
  Update u;
  u.file = file;
  u.mask_file = mask_file;
  u.labels_file = labels_file;
  u.colors = colors;
  u.dOpacity = dOpacity;
  m_mapPendingUpdates[n] = u;
  StartNextUpdate();
}

void QThumbnailWidget::OnThreadFinished()
{
  // finished() is emitted just before the thread ends, so give it a moment.
  // Still running after that: a later run was started, and its own finished()
  // drains the queue.
  if (m_thread->wait(20))
    StartNextUpdate();
}

void QThumbnailWidget::StartNextUpdate()
{
  // Also waits for the initial load, so an update always lands after its item exists.
  if (m_thread->isRunning() || m_mapPendingUpdates.isEmpty())
    return;
  const int n = m_mapPendingUpdates.firstKey();
  const Update u = m_mapPendingUpdates.take(n);
  if (!u.labels_file.isEmpty())
    m_thread->UpdateComponents(n, u.file, u.mask_file, u.labels_file, u.colors, u.dOpacity);
  else if (u.file.isEmpty())
    m_thread->UpdateImage(n, u.image);
  else
    m_thread->UpdateFiles(n, u.file, u.mask_file);
}
