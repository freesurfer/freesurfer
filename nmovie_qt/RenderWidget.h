#ifndef RENDERWIDGET_H
#define RENDERWIDGET_H

#include <QtOpenGL/QGLWidget>
#include <QTimer>
#include <QImage>

class RenderWidget : public QWidget
{
    Q_OBJECT
public:
    explicit RenderWidget(QWidget *parent = 0);

    void paintEvent(QPaintEvent *);
    void resizeEvent(QResizeEvent *);

    void mousePressEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);

    int LoadImages(const QStringList& filenames);


    int GetNumberOfImages()
    {
      return m_images.size();
    }

signals:
    void CurrentImageChanged(const QImage& image, int nIndex );
    void ErrorMessage(const QString& msg);

public slots:
  void OnLoop();
  void OnSwing();
  void OnStop();
  void OnBack();
  void OnForward();
  void SetSpeed(int n);
  void SetAutoResize(bool bAuto);
  void OnTimer();
  void SetCurrentImageIndex(int n);

private:
  QTimer m_timerRender;
  bool   m_bAutoResize;
  int    m_nCurrentImageIndex;
  QList<QImage> m_images;
  QList<QImage> m_resizedImages;
  bool   m_bSwing;
  int    m_nInterval;
  int    m_nY;
  bool   m_bPressed;
};

#endif // RENDERWIDGET_H
