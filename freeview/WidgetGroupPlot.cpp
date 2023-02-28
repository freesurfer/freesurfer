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


#include "WidgetGroupPlot.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "MyUtils.h"
#include "FSGroupDescriptor.h"
#include <QSettings>

#define DISTANCE_THRESHOLD2   36
#define POINT_RADIUS          7

WidgetGroupPlot::WidgetGroupPlot(QWidget *parent) :
  QWidget(parent), m_bAutoScale(true),
  m_fsgd(NULL),
  m_nCurrentVariableIndex(0),
  m_nCurrentDataIndex(-1),
  m_nPlotType(Point)
{
  setFocusPolicy(Qt::StrongFocus);

  QSettings s;
  restoreGeometry(s.value("FsgdWindow/Geometry").toByteArray());
}

WidgetGroupPlot::~WidgetGroupPlot()
{
  QSettings s;
  s.setValue("FsgdWindow/Geometry", saveGeometry());
}

void WidgetGroupPlot::SetFsgdData(FSGroupDescriptor *fsgd)
{
  m_fsgd = fsgd;
  m_dTR = fsgd->m_dXDelta;
  if (m_dTR <= 0)
    m_dTR = 1;
  m_nCurrentVariableIndex = 0;
  update();
}

void WidgetGroupPlot::paintEvent(QPaintEvent *e)
{
  if (!m_fsgd)
  {
    QWidget::paintEvent(e);
    return;
  }

  QPainter p(this);
  QRectF rc_plot = rect();
  int nMargin = 0;
  rc_plot.adjust(nMargin, nMargin, -nMargin, -nMargin);
  rc_plot.adjust(25, 15, -15, -32);

  if (m_fsgd->m_classes.isEmpty())
  {
    p.fillRect(rc_plot, Qt::black);
    return;
  }

  //  if (!m_fsgd->m_title.isEmpty())
  //  {
  //    rc_plot.adjust(0, 5, 0, 0);
  //    p.setPen(Qt::black);
  //    QRectF rc = rect();
  //    rc.setBottom(rc_plot.top());
  //    rc.adjust(0, 5, 0, 0);
  //    p.drawText(rc, Qt::AlignCenter, m_fsgd->m_title);
  //  }

  QFont fnt = font();
  fnt.setPixelSize(11);
  int nTextLen = qMax(QFontMetrics(fnt).horizontalAdvance(QString::number(m_fsgd->m_variables[0].range[0])),
      QFontMetrics(fnt).horizontalAdvance(QString::number(m_fsgd->m_variables[0].range[1])));
  rc_plot.adjust(nTextLen+6, 0, 0, 0);
  p.fillRect(rc_plot.adjusted(-1, -1, 1, 1), Qt::white);

  p.setPen(Qt::black);

  int nStart = 0, nEnd = m_fsgd->m_variables.size()-1;
  if (m_nCurrentVariableIndex >= 0)
    nStart = nEnd = m_nCurrentVariableIndex;
  for (int n = nStart; n <= nEnd ; n++)
  {
    double dMin = m_fsgd->m_dMeasurementRange[0],
        dMax = m_fsgd->m_dMeasurementRange[1],
        dMin_x = m_fsgd->m_variables[n].range[0],
        dMax_x = m_fsgd->m_variables[n].range[1];
    if (m_bAutoScale)
    {
      double val = (dMax-dMin)/20;
      dMin -= val;
      dMax += val;
      val = (dMax_x-dMin_x)/20;
      dMin_x -= val;
      dMax_x += val;
    }
    if (m_fsgd->m_dMeasurementRange[0] >= 0 && dMin < 0)
      dMin = 0;
    if (dMin == dMax)
      dMax += 1;
    if (m_fsgd->m_variables[n].range[0] >= 0 && dMin_x < 0)
      dMin_x = 0;
    QList<double> data_x, data_y;
    for (int i = 0; i < m_fsgd->m_data.size(); i++)
    {
      data_x << m_fsgd->m_data[i].variable_values[n];
      data_y << m_fsgd->m_data[i].measurement;
    }

    QPointF* pts = new QPointF[data_y.size()];
    QString tooltip_strg;
    for (int i = 0; i < data_y.size(); i++)
    {
      pts[i] = QPointF(rc_plot.left() + (data_x[i]-dMin_x)/(dMax_x-dMin_x)*rc_plot.width(),
                       rc_plot.bottom() - (data_y[i]-dMin)/(dMax-dMin)*rc_plot.height());
      QColor c = m_fsgd->m_classes[m_fsgd->m_data[i].class_id].color;
      if (m_nPlotType == Point)
      {
        DrawMarker(&p, pts[i], m_fsgd->m_classes[m_fsgd->m_data[i].class_id].marker.toLower(),
            c, m_nCurrentDataIndex == i);
      }
      // draw current values
      if (m_nCurrentDataIndex == i)
      {
        tooltip_strg = QString("Subject: %1\nClass: %2")
            .arg(m_fsgd->m_data[i].subject_id)
            .arg(m_fsgd->m_classes[m_fsgd->m_data[i].class_id].label);
        for (int m = n; m <=n; m++)
        {
          tooltip_strg += QString("\n%1: %2").arg(m_fsgd->m_variables[m].label)
              .arg(m_fsgd->m_data[i].variable_values[m]);
        }
        tooltip_strg += QString("\nexternal: %1").arg(m_fsgd->m_data[i].measurement);
      }
    }

    // draw cursor
    /*
    p.setPen(QPen(QBrush(Qt::red), 2));
    p.drawLine(pts[m_nCurrentFrame] - QPointF(0, 20),
               pts[m_nCurrentFrame] + QPointF(0, qMin(20., rc_plot.bottom()-pts[m_nCurrentFrame].y())));
               */
    delete[] pts;

    // draw Y metrics
    p.setPen(QPen(Qt::black));
    p.setFont(fnt);
    double nMetricInterval = 40;
    double dMetricStep =  (dMax - dMin) / (rc_plot.height() / nMetricInterval);
    dMetricStep = MyUtils::RoundToGrid( dMetricStep );
    double dMetricPos = (int)(dMin/dMetricStep)*dMetricStep;
    double y = rc_plot.bottom()-(dMetricPos-dMin)/(dMax-dMin)*rc_plot.height();
    while (y > rc_plot.top())
    {
      if (y <= rc_plot.bottom())
      {
        QString strg = QString::number(dMetricPos);
        if (qAbs(dMetricPos) < 1e-8)
          strg = "0";
        p.drawText(QRectF(rect().left(), y-10, rc_plot.left()-5-rect().left(), 20),
                   Qt::AlignVCenter | Qt::AlignRight, strg);
        p.drawLine(QPoint(rc_plot.left()-2, y), QPoint(rc_plot.left(), y));
      }
      dMetricPos += dMetricStep;
      y = rc_plot.bottom()-(dMetricPos-dMin)/(dMax-dMin)*rc_plot.height();
    }

    p.save();
    p.translate(0, rc_plot.top()+rc_plot.height()/2);
    p.rotate(-90);
    p.drawText(QRect(-100, 0, 200, 20), Qt::AlignCenter, m_fsgd->m_measureName);
    p.restore();

    // draw X metrics
    if (true)
    {
      nMetricInterval = 60;
      dMetricStep =  (dMax_x-dMin_x) / (rc_plot.width() / nMetricInterval);
      dMetricStep = MyUtils::RoundToGrid( dMetricStep );
      dMetricPos = (int)(dMin_x/dMetricStep)*dMetricStep;;
      double x = rc_plot.left() + (dMetricPos-dMin_x)/(dMax_x-dMin_x)*rc_plot.width();
      while (x < rc_plot.right())
      {
        if (x > rc_plot.left())
        {
          QString strg = QString::number(dMetricPos);
          p.drawText(QRectF(x-100, rc_plot.bottom()+5, 200, 20),
                     Qt::AlignTop | Qt::AlignHCenter, strg);

          p.setPen(Qt::black);

        }
        dMetricPos += dMetricStep;
        x = rc_plot.left() + (dMetricPos-dMin_x)/(dMax_x-dMin_x)*rc_plot.width();
        p.drawLine(QPoint(x, rc_plot.bottom()), QPoint(x, rc_plot.bottom()+2));
      }
      if (dMin_x < 0 && dMax_x > 0)
      {
        double y = rc_plot.bottom()-(0-dMin)/(dMax-dMin)*rc_plot.height();
        p.drawLine(QPointF(rc_plot.left(), y), QPointF(rc_plot.right(), y));
      }
    }

    p.setPen(Qt::black);
    p.setBrush(Qt::NoBrush);
    p.drawRect(rc_plot);

    QRectF rc = rect().adjusted(0, 0, 0, -3);
    p.drawText(rc, Qt::AlignBottom | Qt::AlignHCenter, m_fsgd->m_variables[m_nCurrentVariableIndex].label);

    // draw current stats
    if (m_nCurrentDataIndex >= 0)
    {
      rc = p.boundingRect(QRectF(0,0,1,1), Qt::AlignLeft|Qt::AlignVCenter, tooltip_strg);
      int bw = 6;
      rc.adjust(-bw, -bw, bw, bw);
      rc.moveTo(pts[m_nCurrentDataIndex] + QPoint(6, 10));
      QPointF pt_c = rc.center();
      qreal w = rc.width()+4, h = rc.height()+4;
      pt_c.rx() = qMin(rc_plot.right()-w/2, qMax(pt_c.rx(), rc_plot.left()+w/2));
      pt_c.ry() = qMin(rc_plot.bottom()-h/2, qMax(pt_c.ry(), rc_plot.top()+h/2));
      rc.moveCenter(pt_c);
      p.setPen(Qt::black);
      p.setBrush(QBrush(QColor(255, 255, 192, 220)));
      p.drawRect(rc);
      rc.adjust(bw, bw, -bw, -bw);
      p.drawText(rc, Qt::AlignLeft | Qt::AlignVCenter, tooltip_strg);
    }
    m_rectPlot = rc_plot;
  }

  p.setRenderHint(QPainter::Antialiasing, false);
}

void WidgetGroupPlot::DrawMarker(QPainter *p, const QPointF &pt, const QString &marker, const QColor &c, bool bHighlight)
{
  DrawMarker(p, pt, marker, c, POINT_RADIUS, bHighlight);
}

void WidgetGroupPlot::DrawMarker(QPainter *p, const QPointF &pt, const QString &marker, const QColor &c, double r, bool bHighlight)
{
  p->save();
  p->setPen(Qt::NoPen);
  p->setBrush(c);
  p->setRenderHint(QPainter::Antialiasing);
  if (marker == "circle")
  {
    p->drawEllipse(pt, r, r);
    r += 3;
  }
  else if (marker == "plus")
  {
    r -= 2;
    p->setPen(QPen(c, 2));
    p->drawLine(pt-QPointF(r, 0), pt+QPointF(r, 0));
    p->drawLine(pt-QPointF(0, r), pt+QPointF(0, r));
    r += 3;
  }
  else if (marker == "dot" || marker == "point")
  {
    p->setRenderHint(QPainter::Antialiasing);
    p->drawEllipse(pt, 3, 3);
    r = 6;
  }
  else if (marker == "triangle")
  {
    p->setRenderHint(QPainter::Antialiasing);
    r += 1;
    double d1 = r*0.866, d2 = r/2;
    QPointF pts[3];
    pts[0] = pt + QPointF(0, -r);
    pts[1] = pt + QPointF(d1, d2);
    pts[2] = pt + QPointF(-d1, d2);
    p->drawPolygon(pts, 3);
    r += 3;
  }
  else if (marker == "square")
  {
    p->drawRect(QRectF(pt-QPointF(r, r), pt+QPointF(r,r)));
    r += 4;
  }
  else if (marker == "diamond")
  {
    p->setRenderHint(QPainter::Antialiasing);
    QPointF pts[4];
    r += 2;
    pts[0] = pt + QPointF(0, -r);
    pts[1] = pt + QPointF(r/2, 0);
    pts[2] = pt + QPointF(0, r);
    pts[3] = pt + QPointF(-r/2, 0);
    p->drawPolygon(pts, 4);
    r += 2;
  }
  else if (marker == "cross")
  {
    p->setPen(QPen(c, 2));
    r -= 2;
    p->drawLine(pt-QPointF(r, r), pt+QPointF(r, r));
    p->drawLine(pt-QPointF(r, -r), pt+QPointF(r, -r));
    r += 5;
  }
  else if (marker == "asterisk")
  {
    r -= 2;
    p->setPen(QPen(c, 2));
    double d = r*0.866;
    p->drawLine(pt-QPointF(d, r/2), pt+QPointF(d, r/2));
    p->drawLine(pt-QPointF(d, -r/2), pt+QPointF(d, -r/2));
    p->drawLine(pt-QPointF(0, r), pt+QPointF(0, r));
    r += 5;
  }
  if (bHighlight)
  {
    p->setPen(QPen(c,1));
    p->setBrush(Qt::NoBrush);
    p->drawEllipse(pt, r, r);
  }
  p->restore();
}

void WidgetGroupPlot::SetAutoScale(bool bAutoScale)
{
  m_bAutoScale = bAutoScale;
  update();
}

void WidgetGroupPlot::SetCurrentVariableIndex(int n)
{
  if (n >= m_fsgd->m_variables.size())
    m_nCurrentVariableIndex = -1;
  else
    m_nCurrentVariableIndex = n;
  update();
}

void WidgetGroupPlot::SetPlotType(int n)
{
  m_nPlotType = n;
  update();
}

void WidgetGroupPlot::mousePressEvent(QMouseEvent *e)
{
  if (m_rectPlot.contains(e->pos()))
  {
    QPointF c_pt = e->pos();
    m_nCurrentDataIndex = -1;
    int n = m_nCurrentVariableIndex;
    double dMin = m_fsgd->m_dMeasurementRange[0],
        dMax = m_fsgd->m_dMeasurementRange[1],
        dMin_x = m_fsgd->m_variables[n].range[0],
        dMax_x = m_fsgd->m_variables[n].range[1];
    if (m_bAutoScale)
    {
      double val = (dMax-dMin)/20;
      dMin -= val;
      dMax += val;
      val = (dMax_x-dMin_x)/20;
      dMin_x -= val;
      dMax_x += val;
    }
    if (m_fsgd->m_dMeasurementRange[0] >= 0 && dMin < 0)
      dMin = 0;
    if (dMin == dMax)
      dMax += 1;
    if (m_fsgd->m_variables[n].range[0] >= 0 && dMin_x < 0)
      dMin_x = 0;
    for (int i = 0; i < m_fsgd->m_data.size(); i++)
    {
      QPointF pt = QPointF(m_rectPlot.left() + (m_fsgd->m_data[i].variable_values[n]-dMin_x)/(dMax_x-dMin_x)*m_rectPlot.width(),
                           m_rectPlot.bottom() - (m_fsgd->m_data[i].measurement-dMin)/(dMax-dMin)*m_rectPlot.height());
      if (m_nPlotType == Point)
      {
        bool bHit = ((pt.x()-c_pt.x())*(pt.x()-c_pt.x()) +
                     (pt.y()-c_pt.y())*(pt.y()-c_pt.y()) < DISTANCE_THRESHOLD2);
        if (bHit)
          m_nCurrentDataIndex = i;
      }
    }
    emit CurrentDataIndexChanged(m_nCurrentDataIndex);
    update();
  }
}

void WidgetGroupPlot::mouseMoveEvent(QMouseEvent *e)
{
  Q_UNUSED(e);
}

void WidgetGroupPlot::leaveEvent(QEvent *e)
{
  Q_UNUSED(e);
}

void WidgetGroupPlot::keyPressEvent(QKeyEvent *e)
{
  Q_UNUSED(e);
}

void WidgetGroupPlot::SetCurrentVertex(int nVertex)
{
  m_fsgd->UpdateData(nVertex);
  update();
}
