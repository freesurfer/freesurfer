#ifndef RENDERWIDGET_H
#define RENDERWIDGET_H

#include <QGLWidget>
#include <QTimer>
#include <QImage>

class RenderWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit RenderWidget(QWidget *parent = 0);

    void paintEvent(QPaintEvent *);
    void resizeEvent(QResizeEvent *);

    int LoadImages(const QStringList& filenames);

    void SetCurrentImageIndex(int n);

signals:
    void CurrentImageChanged(const QImage& image );
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

private:
  QTimer m_timerRender;
  bool   m_bAutoResize;
  int    m_nCurrentImageIndex;
  QList<QImage> m_images;
  QList<QImage> m_resizedImages;
  bool   m_bSwing;
  int    m_nInterval;
};

#endif // RENDERWIDGET_H
