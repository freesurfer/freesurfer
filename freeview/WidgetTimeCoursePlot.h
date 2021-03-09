/**
 * @brief Widget drawing time course plot
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef WIDGETTIMECOURSEPLOT_H
#define WIDGETTIMECOURSEPLOT_H

#include <QWidget>
#include <QList>

struct TimeCourseData
{
public:
  TimeCourseData()
  {
    m_bShow = true;
    m_dXOffset = 0;
    m_dXInterval = 1;
  }
  QList<double>   m_points;
  double          m_dMin;
  double          m_dMax;
  QString         m_strXUnit;
  double          m_dXOffset;
  double          m_dXInterval;
  qint64          m_nId;
  QString         m_strName;

  bool            m_bShow;
  QColor          m_color;
};

class WidgetTimeCoursePlot : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetTimeCoursePlot(QWidget* parent = 0);
  ~WidgetTimeCoursePlot();

  void paintEvent(QPaintEvent * e);

  void AddTimeCourseData(const TimeCourseData& data);
  void SetDataVisible(qint64 nId, bool bVisible);
  void SetDataColor(qint64 nId, const QColor& color);

  void mousePressEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);
  void enterEvent(QEvent* e);
  void leaveEvent(QEvent* e);

  void GetPlotRange(double* range)
  {
    range[0] = m_dMinPlot;
    range[1] = m_dMaxPlot;
  }

  void SetPlotRange(double* range_in);

  void SetShowFrameNumber(bool b)
  {
    m_bShowFrameNumber = b;
    update();
  }

public slots:
  void SetCurrentFrame(int frame);
  void SetAutoScale(bool bAutoScale);
  void ResetPlotRange();
  void Clear();
  void SetDarkMode(bool bDark);

signals:
  void FrameChanged(int frame);
  void PlotRangeChanged();

private:
  QList<TimeCourseData>   m_data;
  double          m_dMin;
  double          m_dMax;
  double          m_dLocalMin;
  double          m_dLocalMax;
  double          m_dMinPlot;
  double          m_dMaxPlot;
  bool            m_bAutoScale;
  int             m_nCurrentFrame;
  int             m_nFrames;
  QRectF          m_rectPlot;
  bool            m_bShowFrameNumber;

  QColor          m_colorBackground;
  QColor          m_colorForeground;
  bool            m_bDarkMode;

  bool            m_bShowCursorInfo;
};

#endif // WIDGETTIMECOURSEPLOT_H
