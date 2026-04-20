#ifndef QTHUMBNAILWIDGET_H
#define QTHUMBNAILWIDGET_H

#include <QListWidget>
#include <QList>
#include <QImage>
#include <QPixmap>

class MakeThumbnailThread;

class QThumbnailWidget : public QListWidget
{
  Q_OBJECT
public:
  QThumbnailWidget(QWidget* parent = NULL);
  virtual ~QThumbnailWidget();

  void LoadImages(const QList<QImage>& images);
  void LoadImageFiles(const QStringList& files, const QStringList& mask_files = QStringList());
  void SetItemImage(int n, const QImage& image);

signals:
  void CurrentIndexChanged(int n);

public slots:
  void SetCurrentIndex(int n);
  void SetItemChecked(int n, bool bChecked = true);
  void OnPixmapReady(const QPixmap& pix, int n);
  void Clear();

protected slots:
  void OnItemSelectionChanged();

private:
  QList<QPixmap>  m_listPixmap;
  MakeThumbnailThread* m_thread;
};

#endif // QTHUMBNAILWIDGET_H
