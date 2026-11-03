#ifndef QTHUMBNAILWIDGET_H
#define QTHUMBNAILWIDGET_H

#include <QListWidget>
#include <QList>
#include <QImage>
#include <QPixmap>
#include <QMap>
#include <QSet>
#include <QColor>

class MakeThumbnailThread;

class QThumbnailWidget : public QListWidget
{
  Q_OBJECT
public:
  QThumbnailWidget(QWidget* parent = NULL);
  virtual ~QThumbnailWidget();

  // Extra space per axis around a thumbnail's icon (3px border + 3px icon-mode
  // margin, per side). Used by ImageGridView's cell maths.
  static int ItemPadding();

  // 'checked_indices': which 'files' get the green "processed" outline as each
  // thumbnail loads (e.g. photos already corrected when Correction is
  // re-entered).
  void LoadImageFiles(const QStringList& files, const QStringList& mask_files = QStringList(),
                       const QList<int>& checked_indices = QList<int>());
  void SetItemImage(int n, const QImage& image);
  // Rebuilds item n from a photo and its mask on disk, off the GUI thread.
  void SetItemFiles(int n, const QString& file, const QString& mask_file);
  // Rebuilds item n as the Connected components view shows the photo (see
  // MakeThumbnailThread::UpdateComponents()), off the GUI thread.
  void SetItemComponents(int n, const QString& file, const QString& mask_file, const QString& labels_file,
                         const QList<QColor>& colors, double dOpacity);

signals:
  void CurrentIndexChanged(int n);

public slots:
  void SetCurrentIndex(int n);
  void SetItemChecked(int n, bool bChecked = true);
  void OnPixmapReady(const QPixmap& pix, int n, int nGeneration);
  void Clear();

protected slots:
  void OnItemSelectionChanged();
  void OnThreadFinished();

private:
  struct Update { QImage image; QString file, mask_file, labels_file; QList<QColor> colors; double dOpacity = 1; };
  void StartNextUpdate();

  QList<QPixmap>  m_listPixmap;
  // Item updates wait here while the thread is busy (QThread::start() does
  // nothing on a running thread); the latest one per item wins.
  QMap<int, Update> m_mapPendingUpdates;
  MakeThumbnailThread* m_thread;
  QList<int>      m_listPendingChecked;
  // Checked items, re-marked whenever their thumbnail is rebuilt.
  QSet<int>       m_setChecked;
  // Bumped by Clear(): thumbnails already queued for the old list are dropped.
  int             m_nGeneration = 0;
};

#endif // QTHUMBNAILWIDGET_H
