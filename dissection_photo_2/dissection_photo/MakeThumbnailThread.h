#ifndef MAKETHUMBNAILTHREAD_H
#define MAKETHUMBNAILTHREAD_H

#include <QPixmap>
#include <QThread>
#include <QMutex>

class MakeThumbnailThread : public QThread
{
  Q_OBJECT
public:
  MakeThumbnailThread(QObject* parent = NULL);

signals:
  void PixmapReady(const QPixmap& pix, int n = -1);

public slots:
  void LoadImages(const QList<QImage>& images)
  {
    m_images = images;
    m_files.clear();
    m_maskFiles.clear();
    m_nUpdateIndex = -1;
    start();
  }

  void LoadImageFiles(const QStringList& files, const QStringList& mask_files = QStringList())
  {
    m_files = files;
    m_maskFiles = mask_files;
    m_images.clear();
    m_nUpdateIndex = -1;
    start();
  }

  void UpdateImage(int n, const QImage& img);

  void SetPixmapSize(const QSize& sz, int devicePixelRatio)
  {
    m_size = sz;
    m_ratio = devicePixelRatio;
  }

  void Abort();

protected:
  void run();
  QPixmap BuildPixmap(int nIndex, const QImage& img, const QImage& mask = QImage());

  QMutex          m_mutex;
  QList<QImage>   m_images;
  QStringList     m_files;
  QStringList     m_maskFiles;
  QSize         m_size;
  int           m_ratio;
  bool          m_bAbort;
  int           m_nUpdateIndex;
};

#endif // MAKETHUMBNAILTHREAD_H
