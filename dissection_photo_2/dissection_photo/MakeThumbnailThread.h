#ifndef MAKETHUMBNAILTHREAD_H
#define MAKETHUMBNAILTHREAD_H

#include <QPixmap>
#include <QColor>
#include <QList>
#include <QThread>
#include <QMutex>

class MakeThumbnailThread : public QThread
{
  Q_OBJECT
public:
  MakeThumbnailThread(QObject* parent = NULL);

signals:
  // nGeneration: the SetGeneration() value the work was started under, so the
  // receiver can drop results of a list it has since replaced.
  void PixmapReady(const QPixmap& pix, int n, int nGeneration);

public slots:
  void LoadImageFiles(const QStringList& files, const QStringList& mask_files = QStringList())
  {
    m_files = files;
    m_maskFiles = mask_files;
    m_images.clear();
    m_nUpdateIndex = -1;
    start();
  }

  // One thumbnail, from an image or from a photo and its mask on disk. Only
  // while the thread is idle: QThumbnailWidget queues them.
  void UpdateImage(int n, const QImage& img);
  void UpdateFiles(int n, const QString& file, const QString& mask_file);
  // A Connected components thumbnail: the photo inside its mask, with the saved
  // components (a label .npy) in colour on top, as the view shows them.
  void UpdateComponents(int n, const QString& file, const QString& mask_file, const QString& labels_file,
                        const QList<QColor>& colors, double dOpacity);

  void SetPixmapSize(const QSize& sz, int devicePixelRatio)
  {
    m_size = sz;
    m_ratio = devicePixelRatio;
  }

  void Abort();
  // Set only while the thread is not running.
  void SetGeneration(int n) { m_nGeneration = n; }

protected:
  void run();
  QPixmap BuildPixmap(int nIndex, const QImage& img, const QImage& mask = QImage());
  QPixmap BuildPixmapFromFiles(int nIndex, const QString& file, const QString& mask_file);
  QPixmap BuildComponentsPixmap(int nIndex);

  QMutex          m_mutex;
  QList<QImage>   m_images;
  QStringList     m_files;
  QStringList     m_maskFiles;
  QSize         m_size;
  int           m_ratio;
  bool          m_bAbort;
  int           m_nUpdateIndex;
  int           m_nGeneration = 0;
  QString       m_strUpdateFile;
  QString       m_strUpdateMask;
  QString       m_strUpdateLabels;
  QList<QColor> m_listColors;
  double        m_dOpacity = 1;
};

#endif // MAKETHUMBNAILTHREAD_H
