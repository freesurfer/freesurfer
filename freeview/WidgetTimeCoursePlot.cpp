/**
 * @file  WidgetTimeCoursePlot.cpp
 * @brief Widget drawing time course plot
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2017/02/01 15:28:54 $
 *    $Revision: 1.9 $
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
  QWidget(parent), m_bAutoScale(true), m_nCurrentFrame(0), m_dMinPlot(0), m_dMaxPlot(1),
  m_dXInterval(1), m_dXOffset(0), m_bShowFrameNumber(false)
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

  if (m_dMin != min_val || m_dMax != max_val )
  {
    m_dMinPlot = min_val;
    m_dMaxPlot = max_val;
    emit PlotRangeChanged();
  }
  m_dMin = min_val;
  m_dMax = max_val;
  if (m_dMax <= m_dMin)
    m_dMax = m_dMin+1;
  if (m_nCurrentFrame >= data.size())
    m_nCurrentFrame = 0;
  update();
}

void WidgetTimeCoursePlot::SetSecondData(const QList<double> &data)
{
  m_secondData = data;
  update();
}

void WidgetTimeCoursePlot::paintEvent(QPaintEvent *e)
{
  Q_UNUSED(e);
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

  double dMin = m_dMinPlot, dMax = m_dMaxPlot;
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
    if (!m_secondData.isEmpty())
    {
      for (int i = 0; i < m_secondData.size(); i++)
      {
        if (dMin > m_secondData[i])
          dMin = m_secondData[i];
        else if (dMax < m_secondData[i])
          dMax = m_secondData[i];
      }
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
  QPointF* pts2 = NULL;
  if (m_secondData.size() == m_data.size())
    pts2 = new QPointF[m_data.size()];
  for (int i = 0; i < m_data.size(); i++)
  {
    pts[i] = QPointF(rc_plot.left() + dSpacing*i,
                     rc_plot.bottom() - (m_data[i]-dMin)/(dMax-dMin)*rc_plot.height());
    if (m_secondData.size() >= m_data.size())
      pts2[i] = QPointF(rc_plot.left() + dSpacing*i,
                       rc_plot.bottom() - (m_secondData[i]-dMin)/(dMax-dMin)*rc_plot.height());
  }
  p.setPen(QPen(QBrush(Qt::yellow), 2));
  p.setClipRect(rc_plot);
  p.drawPolyline(pts, m_data.size());
  if (pts2)
  {
    p.setPen(QPen(QBrush(Qt::cyan), 2));
    p.drawPolyline(pts2, m_data.size());
  }

  // draw cursor
  p.setPen(QPen(QBrush(Qt::red), 2));
  p.drawLine(pts[m_nCurrentFrame].x(), rc_plot.top(), pts[m_nCurrentFrame].x(), rc_plot.bottom());

  p.setClipping(false);
  delete[] pts;
  delete[] pts2;

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
      p.drawText(QRectF(rect().left(), y-10, rc_plot.left()-rect().left()-6, 20),
                 Qt::AlignVCenter | Qt::AlignRight, strg);
      p.drawLine(QPointF(rc_plot.left(), y), QPointF(rc_plot.left()-2, y));
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
  dMetricStep =  (m_data.size()-1) / (rc_plot.width() / nMetricInterval);
  dMetricStep = MyUtils::RoundToGrid( dMetricStep );
  dMetricPos = 0;
  double dScale = 1;
  int nPrecise = 0;
  double x = rc_plot.left();
  QString strXUnit = m_strXUnit;
  if (!m_bShowFrameNumber)
  {
    double np = log10(qAbs(m_dXInterval*dMetricStep));
    if (np > 3 || np < -3)
    {
      np = ((int)np);
      dScale = pow(10, -np);
      strXUnit = QString("10^%1 %2").arg(np).arg(m_strXUnit);
      nPrecise = 2;
    }
    else if (np < 2 && qAbs(m_dXInterval*dMetricStep) < 10)
      nPrecise = 2;

    dMetricStep = (m_data.size()-1) * qAbs(m_dXInterval) / rc_plot.width() * nMetricInterval;
    dMetricStep = MyUtils::RoundToGrid(dMetricStep)/qAbs(m_dXInterval);
    double dval = m_dXOffset/qAbs(m_dXInterval*dMetricStep) + 1;
    if (m_dXInterval > 0)
      dval = ((int)dval) + 1 - dval;
    else
      dval = dval - ((int)dval);
    dMetricPos = dval*dMetricStep - dMetricStep;
    while (dMetricPos < -dMetricStep/10)
      dMetricPos += dMetricStep;
    x = rc_plot.left() + dMetricPos*rc_plot.width()/(m_data.size()-1);
  }
  while (x <= rc_plot.right()+dMetricStep/10)
  {
    QString strg;
    if (m_bShowFrameNumber)
      strg = QString::number(dMetricPos);
    else
      strg = QString::number((m_dXOffset+m_dXInterval*dMetricPos)*dScale, 'f', nPrecise);
    p.drawText(QRectF(x-100, rc_plot.bottom()+5, 200, 20),
               Qt::AlignTop | Qt::AlignHCenter, strg);
    if (x-1 >= rc_plot.left() && x-1 <= rc_plot.right())
      p.drawLine(QPointF(x-1, rc_plot.bottom()), QPointF(x-1, rc_plot.bottom()+2));

    dMetricPos += dMetricStep;
    x = rc_plot.left() + dMetricPos*rc_plot.width()/(m_data.size()-1);
  }

  QRectF rc = rect();
  // draw current x value
  QString x_strg = QString("Frame #%1").arg(m_nCurrentFrame);
  if (!strXUnit.isEmpty())
    x_strg += QString(" / %1 (%2)").arg((m_dXOffset+m_nCurrentFrame*m_dXInterval)*dScale).arg(strXUnit);
  p.drawText(rc, Qt::AlignBottom | Qt::AlignHCenter, x_strg);

  // draw current y value
  QString strg = QString("Signal intensity: %1 %2")
      .arg(m_data[m_nCurrentFrame])
      .arg(!m_secondData.isEmpty() ? QString("/ %1").arg(m_secondData[m_nCurrentFrame]) : "");
  rc = rect().adjusted(0, 5, 0, 0);
  p.drawText(rc, Qt::AlignHCenter | Qt::AlignTop, strg);
  m_rectPlot = rc_plot;
  if (m_dMinPlot != dMin || m_dMaxPlot != dMax)
  {
    m_dMinPlot = dMin;
    m_dMaxPlot = dMax;
    emit PlotRangeChanged();
  }
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
  if (e->button() == Qt::LeftButton && m_rectPlot.contains(e->pos()))
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
  if (e->buttons() & Qt::LeftButton && m_rectPlot.contains(e->pos()))
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
  QWidget::keyPressEvent(e);
}

void WidgetTimeCoursePlot::SetPlotRange(double *range_in)
{
  m_dMinPlot = range_in[0];
  m_dMaxPlot = range_in[1];
  update();
}

void WidgetTimeCoursePlot::ResetPlotRange()
{
  m_dMinPlot = m_dMin;
  m_dMaxPlot = m_dMax;
  update();
  emit PlotRangeChanged();
}
