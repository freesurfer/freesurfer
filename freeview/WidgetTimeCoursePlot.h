/**
 * @file  WidgetTimeCoursePlot.h
 * @brief Widget drawing time course plot
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/08/29 15:24:59 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

class WidgetTimeCoursePlot : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetTimeCoursePlot(QWidget* parent = 0);
  ~WidgetTimeCoursePlot();

  void paintEvent(QPaintEvent * e);

  void SetTimeCourseData(const QList<double>& data, double min_val, double max_val,
                         double t_interval = 1000);
  void SetSecondData(const QList<double>& data);

  void mousePressEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);

  void GetPlotRange(double* range)
  {
    range[0] = m_dMinPlot;
    range[1] = m_dMaxPlot;
  }

  void SetPlotRange(double* range_in);

  void SetXUnitInfo(double interval, double offset, const QString& unit_str = "")
  {
    m_dXInterval = interval;
    m_dXOffset = offset;
    m_strXUnit = unit_str;
    update();
  }

  void SetShowFrameNumber(bool b)
  {
    m_bShowFrameNumber = b;
    update();
  }

public slots:
  void SetCurrentFrame(int frame);
  void SetAutoScale(bool bAutoScale);
  void ResetPlotRange();

signals:
  void FrameChanged(int frame);
  void PlotRangeChanged();

private:
  QList<double>   m_data;
  QList<double>   m_secondData;
  double          m_dTR;
  double          m_dMin;
  double          m_dMax;
  double          m_dMinPlot;
  double          m_dMaxPlot;
  bool            m_bAutoScale;
  int             m_nCurrentFrame;
  QRectF          m_rectPlot;

  QString         m_strXUnit;
  double          m_dXOffset;
  double          m_dXInterval;
  bool            m_bShowFrameNumber;
};

#endif // WIDGETTIMECOURSEPLOT_H
