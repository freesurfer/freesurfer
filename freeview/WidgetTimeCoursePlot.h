#ifndef WIDGETTIMECOURSEPLOT_H
#define WIDGETTIMECOURSEPLOT_H

#include <QWidget>
#include <QList>

class WidgetTimeCoursePlot : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetTimeCoursePlot(QWidget* parent = 0);
  ~WidgetTimeCoursePlot();

  void paintEvent(QPaintEvent * e);

  void SetTimeCourseData(const QList<double>& data, double min_val, double max_val,
                         double t_interval = -1);

  void mousePressEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);

public slots:
  void SetCurrentFrame(int frame);
  void SetAutoScale(bool bAutoScale);

signals:
  void FrameChanged(int frame);

private:
  QList<double>   m_data;
  double          m_dTR;
  double          m_dMin;
  double          m_dMax;
  bool            m_bAutoScale;
  int             m_nCurrentFrame;
  QRectF          m_rectPlot;
};

#endif // WIDGETTIMECOURSEPLOT_H
