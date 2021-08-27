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


#include "WidgetTimeCoursePlot.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "MyUtils.h"

WidgetTimeCoursePlot::WidgetTimeCoursePlot(QWidget *parent) :
  QWidget(parent), m_bAutoScale(true), m_nCurrentFrame(0), m_dMinPlot(0), m_dMaxPlot(1),
  m_bShowFrameNumber(false), m_bShowCursorInfo(false)
{
  setFocusPolicy(Qt::StrongFocus);
  SetDarkMode(true);
}

WidgetTimeCoursePlot::~WidgetTimeCoursePlot()
{
}

void WidgetTimeCoursePlot::SetDarkMode(bool bDark)
{
  m_bDarkMode = bDark;
  m_colorBackground = bDark?Qt::black:Qt::white;
  m_colorForeground = bDark?QColor(255,255,255):QColor(30,30,30);
}

void WidgetTimeCoursePlot::Clear()
{
  m_data.clear();
  update();
}

void WidgetTimeCoursePlot::AddTimeCourseData(const TimeCourseData &data)
{
  m_data << data;
  if (m_data.size() == 1) // first one
  {
    m_nFrames = 0;
    m_dMin = m_dMax = m_dLocalMin = m_dLocalMax = data.m_points.first();
  }

  if (data.m_points.size() > m_nFrames)
    m_nFrames = data.m_points.size();
  if (m_dMin > data.m_dMin)
    m_dMin = data.m_dMin;
  if (m_dMax < data.m_dMax)
    m_dMax = data.m_dMax;

  for (int i = 0; i < data.m_points.size(); i++)
  {
    if (m_dLocalMin > data.m_points[i])
      m_dLocalMin = data.m_points[i];
    else if (m_dLocalMax < data.m_points[i])
      m_dLocalMax = data.m_points[i];
  }

  if (m_dMax <= m_dMin)
    m_dMax = m_dMin+1;

  if (m_dLocalMax <= m_dLocalMin)
    m_dLocalMax = m_dLocalMin+1;

  if (m_nCurrentFrame >= data.m_points.size())
    m_nCurrentFrame = 0;

  if (m_dMinPlot == 0 && m_dMaxPlot == 1)
  {
    m_dMinPlot = m_dMin;
    m_dMaxPlot = m_dMax;
  }

  update();
}

void WidgetTimeCoursePlot::SetDataVisible(qint64 nId, bool bShow)
{
  for (int i = 0; i < m_data.size(); i++)
  {
    if (m_data[i].m_nId == nId)
    {
      m_data[i].m_bShow = bShow;
      update();
      break;
    }
  }
}

void WidgetTimeCoursePlot::SetDataColor(qint64 nId, const QColor& color)
{
  for (int i = 0; i < m_data.size(); i++)
  {
    if (m_data[i].m_nId == nId)
    {
      m_data[i].m_color = color;
      update();
      break;
    }
  }
}


void WidgetTimeCoursePlot::paintEvent(QPaintEvent *e)
{
  Q_UNUSED(e);
  QPainter p(this);
  QRectF rc_plot = rect();
  p.fillRect(rect(), m_colorBackground);
  int nMargin = 10;
  rc_plot.adjust(nMargin, nMargin, -nMargin, -nMargin);
  rc_plot.adjust(15, 18, -25, -30);

  if (m_data.isEmpty())
  {
    return;
  }

  QFont fnt = font();
  fnt.setPixelSize(11);
  int nTextLen = qMax(QFontMetrics(fnt).width(QString::number(m_dMax)),
                      QFontMetrics(fnt).width(QString::number(m_dMin)));
  rc_plot.adjust(nTextLen+6, 0, 0, 0);
  p.fillRect(rc_plot.adjusted(-1, -1, 1, 1), m_colorBackground);

  double dMin = m_dMinPlot, dMax = m_dMaxPlot;
  if (m_bAutoScale)
  {
    dMin = m_dLocalMin;
    dMax = m_dLocalMax;
    dMax += (dMax-dMin)/4;
    double old_min = dMin;
    dMin -= (dMax-dMin)/4;
    if (dMin < 0 && old_min >= 0)
      dMin = 0;
  }
  if (dMin == dMax)
    dMax += 1;

  p.setRenderHint(QPainter::Antialiasing);
  // draw rect
  p.setPen(QPen(m_colorForeground, 1));
  p.setBrush(Qt::NoBrush);
  p.drawRect(rc_plot);

  // draw plots
  double dSpacing = rc_plot.width() / (m_nFrames-1);
  p.save();
  for (int n = 0; n < m_data.size(); n++)
  {
    TimeCourseData& td = m_data[n];
    if (td.m_bShow)
    {
      QPointF* pts = new QPointF[td.m_points.size()];
      for (int i = 0; i < td.m_points.size(); i++)
      {
        pts[i] = QPointF(rc_plot.left() + dSpacing*i,
                         rc_plot.bottom() - (td.m_points[i]-dMin)/(dMax-dMin)*rc_plot.height());
      }
      p.setPen(QPen(QBrush(td.m_color), 2));
      p.setClipRect(rc_plot);
      p.drawPolyline(pts, td.m_points.size());
      delete[] pts;
    }
  }
  p.restore();

  // draw cursor
  p.setPen(QPen(QColor(255,255,255), 1));
  int cx = rc_plot.left()+dSpacing*m_nCurrentFrame;
  p.drawLine(cx, rc_plot.top(), cx, rc_plot.bottom());
  p.setClipping(false);

  // draw Y metrics
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
  dMetricStep =  (m_nFrames-1) / (rc_plot.width() / nMetricInterval);
  dMetricStep = MyUtils::RoundToGrid( dMetricStep );
  dMetricPos = 0;
  double dScale = 1;
  int nPrecise = 0;
  double x = rc_plot.left();
  TimeCourseData& td = m_data[0];
  QString strXUnit = td.m_strXUnit;
  if (!m_bShowFrameNumber)
  {
    double np = log10(qAbs(td.m_dXInterval*dMetricStep));
    if (np > 3 || np < -3)
    {
      np = ((int)np);
      dScale = pow(10, -np);
      strXUnit = QString("10^%1 %2").arg(np).arg(td.m_strXUnit);
      nPrecise = 2;
    }
    else if (np < 2 && qAbs(td.m_dXInterval*dMetricStep) < 10)
      nPrecise = 2;

    dMetricStep = (m_nFrames-1) * qAbs(td.m_dXInterval) / rc_plot.width() * nMetricInterval;
    dMetricStep = MyUtils::RoundToGrid(dMetricStep)/qAbs(td.m_dXInterval);
    double dval = td.m_dXOffset/qAbs(td.m_dXInterval*dMetricStep) + 1;
    if (td.m_dXInterval > 0)
      dval = ((int)dval) + 1 - dval;
    else
      dval = dval - ((int)dval);
    dMetricPos = dval*dMetricStep - dMetricStep;
    while (dMetricPos < -dMetricStep/10)
      dMetricPos += dMetricStep;
    x = rc_plot.left() + dMetricPos*rc_plot.width()/(m_nFrames-1);
  }
  while (x <= rc_plot.right()+dMetricStep/10)
  {
    QString strg;
    if (m_bShowFrameNumber)
      strg = QString::number(dMetricPos);
    else
      strg = QString::number((td.m_dXOffset+td.m_dXInterval*dMetricPos)*dScale, 'f', nPrecise);
    p.drawText(QRectF(x-100, rc_plot.bottom()+5, 200, 20),
               Qt::AlignTop | Qt::AlignHCenter, strg);
    if (x-1 >= rc_plot.left() && x-1 <= rc_plot.right())
      p.drawLine(QPointF(x-1, rc_plot.bottom()), QPointF(x-1, rc_plot.bottom()+2));

    dMetricPos += dMetricStep;
    x = rc_plot.left() + dMetricPos*rc_plot.width()/(m_nFrames-1);
  }

  QRectF rc = rect();
  // draw current x value
  QString x_strg = QString("Frame #%1").arg(m_nCurrentFrame);
  if (!strXUnit.isEmpty())
    x_strg += QString(" / %1 (%2)").arg((td.m_dXOffset+m_nCurrentFrame*td.m_dXInterval)*dScale).arg(strXUnit);
  p.drawText(rc, Qt::AlignBottom | Qt::AlignHCenter, x_strg);

  // draw current y values
  if (m_bShowCursorInfo)
  {
    int nMaxNameWidth = 0, nMaxValueWidth = 0;
    QFontMetrics fmt(p.font());
    int nVisibleLines = 0;
    for (int n = 0; n < m_data.size(); n++)
    {
      TimeCourseData& td = m_data[n];
      if (td.m_bShow)
      {
        int nLen = fmt.width(td.m_strName+" :");
        if (nLen > nMaxNameWidth)
          nMaxNameWidth = nLen;
        if (m_nCurrentFrame < td.m_points.size())
        {
          nLen = fmt.width(QString::number(td.m_points[m_nCurrentFrame]));
          if (nLen > nMaxValueWidth)
            nMaxValueWidth = nLen;
        }
        nVisibleLines++;
      }
    }
    if (nMaxValueWidth > 0)
    {
      int nMarginH = 17, nMarginV = 15;
      int nSpacingH = 9, nSpacingV = 5;
      QRectF rc_frame(0, 0, nMarginH*2 + nSpacingH + nMaxValueWidth + nMaxNameWidth,
               nMarginV*2 + (nVisibleLines-1)*nSpacingV + nVisibleLines*fmt.height());
      rc_frame.moveTopRight(rc_plot.topRight() + QPointF(-15, 15));
      p.setPen(QPen(m_colorForeground, 1));
      p.setBrush(m_colorBackground);
      p.drawRect(rc_frame);
      int offset_y = nMarginV;
      for (int n = 0; n < m_data.size(); n++)
      {
        TimeCourseData& td = m_data[n];
        if (td.m_bShow)
        {
          p.setPen(QPen(td.m_color));
          QRectF rc(0, 0, nMaxNameWidth, fmt.height());
          rc.moveTopLeft(rc_frame.topLeft() + QPointF(nMarginH, offset_y));
          p.drawText(rc, Qt::AlignRight, td.m_strName+" :");
          if (m_nCurrentFrame < td.m_points.size())
          {
            rc.setWidth(nMaxValueWidth);
            rc.moveRight(rc_frame.right()-nMarginH-3);
            p.setPen(m_colorForeground);
            p.drawText(rc, Qt::AlignLeft, QString::number(td.m_points[m_nCurrentFrame]));
          }
          offset_y += fmt.height() + nSpacingV;
        }
      }
    }
  }

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
    double dSpacing = m_rectPlot.width() / (m_nFrames-1);
    int n = (int)((e->x() - m_rectPlot.left() ) / dSpacing + 0.5);
    if (n >= 0 && n < m_nFrames)
    {
      m_nCurrentFrame = n;
      update();
      emit FrameChanged(n);
    }
  }
}

void WidgetTimeCoursePlot::mouseMoveEvent(QMouseEvent *e)
{
  m_bShowCursorInfo = true;
  if (e->buttons() & Qt::LeftButton && m_rectPlot.contains(e->pos()))
  {
    double dSpacing = m_rectPlot.width() / (m_nFrames-1);
    int n = (int)((e->x() - m_rectPlot.left() ) / dSpacing + 0.5);
    if (n >= 0 && n < m_nFrames)
    {
      m_nCurrentFrame = n;
      update();
      emit FrameChanged(n);
    }
  }
}

void WidgetTimeCoursePlot::enterEvent(QEvent *e)
{
  m_bShowCursorInfo = true;
  update();
}

void WidgetTimeCoursePlot::leaveEvent(QEvent *e)
{
  m_bShowCursorInfo = false;
  update();
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
    if (m_nCurrentFrame < m_nFrames-1)
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
