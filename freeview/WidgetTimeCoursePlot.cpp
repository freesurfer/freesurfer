/**
 * @file  WidgetTimeCoursePlot.cpp
 * @brief Widget drawing time course plot
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.8 $
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


#include "WidgetTimeCoursePlot.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "MyUtils.h"

WidgetTimeCoursePlot::WidgetTimeCoursePlot(QWidget *parent) :
  QWidget(parent), m_bAutoScale(true), m_nCurrentFrame(0)
{
  setFocusPolicy(Qt::StrongFocus);
}

WidgetTimeCoursePlot::~WidgetTimeCoursePlot()
{
}

void WidgetTimeCoursePlot::SetTimeCourseData(const QList<double> &data,
                                             double min_val, double max_val,
                                             double t_interval)
{
  m_data = data;
  m_dTR = t_interval;
  if (m_dTR <= 0)
    m_dTR = 1000;
  m_dMin = min_val;
  m_dMax = max_val;
  if (m_dMax <= m_dMin)
    m_dMax = m_dMin+1;
  if (m_nCurrentFrame >= data.size())
    m_nCurrentFrame = 0;
  update();
}

void WidgetTimeCoursePlot::paintEvent(QPaintEvent *e)
{
  QPainter p(this);
  QRectF rc_plot = rect();
  int nMargin = 10;
  rc_plot.adjust(nMargin, nMargin, -nMargin, -nMargin);
  rc_plot.adjust(15, 15, -25, -27);

  if (m_data.isEmpty())
  {
    p.fillRect(rc_plot, Qt::black);
    return;
  }

  QFont fnt = font();
  fnt.setPixelSize(11);
  int nTextLen = qMax(QFontMetrics(fnt).width(QString::number(m_dMax)),
                      QFontMetrics(fnt).width(QString::number(m_dMin)));
  rc_plot.adjust(nTextLen+6, 0, 0, 0);
  p.fillRect(rc_plot.adjusted(-1, -1, 1, 1), Qt::black);

  double dMin = m_dMin, dMax = m_dMax;
  if (m_bAutoScale)
  {
    dMin = dMax = m_data[0];
    for (int i = 1; i < m_data.size(); i++)
    {
      if (dMin > m_data[i])
        dMin = m_data[i];
      else if (dMax < m_data[i])
        dMax = m_data[i];
    }
    dMax += (dMax-dMin)/4;
    double old_min = dMin;
    dMin -= (dMax-dMin)/4;
    if (dMin < 0 && old_min >= 0)
      dMin = 0;
  }
  if (dMin == dMax)
    dMax += 1;
  double dSpacing = rc_plot.width() / (m_data.size()-1);
  p.setRenderHint(QPainter::Antialiasing);
  QPointF* pts = new QPointF[m_data.size()];
  for (int i = 0; i < m_data.size(); i++)
  {
    pts[i] = QPointF(rc_plot.left() + dSpacing*i,
                     rc_plot.bottom() - (m_data[i]-dMin)/(dMax-dMin)*rc_plot.height());
  }
  p.setPen(QPen(QBrush(Qt::yellow), 2));
  p.drawPolyline(pts, m_data.size());

  // draw cursor
  p.setPen(QPen(QBrush(Qt::red), 2));
  p.drawLine(pts[m_nCurrentFrame] - QPointF(0, 20),
             pts[m_nCurrentFrame] + QPointF(0, qMin(20., rc_plot.bottom()-pts[m_nCurrentFrame].y())));
  delete[] pts;

  // draw Y metrics
  p.setPen(QPen(Qt::black));
  p.setFont(fnt);
  double nMetricInterval = 30;
  double dMetricStep =  (dMax - dMin) / (rc_plot.height() / nMetricInterval);
  dMetricStep = MyUtils::RoundToGrid( dMetricStep );
  double dMetricPos = (int)(dMin/dMetricStep)*dMetricStep;
  double y = rc_plot.bottom()-(dMetricPos-dMin)/(dMax-dMin)*rc_plot.height();
  while (y > rc_plot.top())
  {
    if (y <= rc_plot.bottom())
    {
      QString strg = QString::number(dMetricPos);
      p.drawText(QRectF(rect().left(), y-10, rc_plot.left()-5-rect().left(), 20),
                 Qt::AlignVCenter | Qt::AlignRight, strg);
    }
    dMetricPos += dMetricStep;
    y = rc_plot.bottom()-(dMetricPos-dMin)/(dMax-dMin)*rc_plot.height();
  }

  p.save();
  p.translate(0, rc_plot.top()+rc_plot.height()/2);
  p.rotate(-90);
  p.drawText(QRect(-100, 0, 200, 20), Qt::AlignCenter, "Signal Intensity");
  p.restore();

  // draw X metrics
  nMetricInterval = 50;
  double dTR = 1; // m_dTR;
  dMetricStep =  (m_data.size()-1)*dTR / (rc_plot.width() / nMetricInterval);
  dMetricStep = MyUtils::RoundToGrid( dMetricStep );
  dMetricPos = 0;
  double x = rc_plot.left();
  while (x < rc_plot.right())
  {
    QString strg = QString::number(dMetricPos);
    p.drawText(QRectF(x-100, rc_plot.bottom()+5, 200, 20),
               Qt::AlignTop | Qt::AlignHCenter, strg);

    dMetricPos += dMetricStep;
    x = rc_plot.left() + dMetricPos/((m_data.size()-1)*dTR)*rc_plot.width();
  }

  QRectF rc = rect().adjusted(0, 0, 0, -3);
  p.drawText(rc, Qt::AlignBottom | Qt::AlignHCenter, "Frame ");

  // draw current stats
  QString strg = QString("Signal intensity:%1   Frame: %3")
                 .arg(m_data[m_nCurrentFrame])
              //   .arg(m_nCurrentFrame*m_dTR/1000)
                 .arg(m_nCurrentFrame);
  rc = rect().adjusted(0, 5, 0, 0);
  p.drawText(rc, Qt::AlignHCenter | Qt::AlignTop, strg);
  m_rectPlot = rc_plot;
}

void WidgetTimeCoursePlot::SetCurrentFrame(int frame)
{
  m_nCurrentFrame = frame;
  update();
}

void WidgetTimeCoursePlot::SetAutoScale(bool bAutoScale)
{
  m_bAutoScale = bAutoScale;
  update();
}

void WidgetTimeCoursePlot::mousePressEvent(QMouseEvent *e)
{
  if (e->button() == Qt::LeftButton && m_rectPlot.contains(e->posF()))
  {
    double dSpacing = m_rectPlot.width() / (m_data.size()-1);
    int n = (int)((e->x() - m_rectPlot.left() ) / dSpacing + 0.5);
    if (n >= 0 && n < m_data.size())
    {
      m_nCurrentFrame = n;
      update();
      emit FrameChanged(n);
    }
  }
}

void WidgetTimeCoursePlot::mouseMoveEvent(QMouseEvent *e)
{
  if (e->buttons() & Qt::LeftButton && m_rectPlot.contains(e->posF()))
  {
    double dSpacing = m_rectPlot.width() / (m_data.size()-1);
    int n = (int)((e->x() - m_rectPlot.left() ) / dSpacing + 0.5);
    if (n >= 0 && n < m_data.size())
    {
      m_nCurrentFrame = n;
      update();
      emit FrameChanged(n);
    }
  }
}

void WidgetTimeCoursePlot::keyPressEvent(QKeyEvent *e)
{
  if (e->key() == Qt::Key_Left)
  {
    if (m_nCurrentFrame > 0)
    {
      m_nCurrentFrame--;
      update();
      emit FrameChanged(m_nCurrentFrame);
    }
  }
  else if (e->key() == Qt::Key_Right)
  {
    if (m_nCurrentFrame < m_data.size()-1)
    {
      m_nCurrentFrame++;
      update();
      emit FrameChanged(m_nCurrentFrame);
    }
  }
}
