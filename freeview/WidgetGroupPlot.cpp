/**
 * @file  WidgetGroupPlot.cpp
 * @brief Widget drawing time course plot
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/04/25 00:04:02 $
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


#include "WidgetGroupPlot.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDebug>
#include "MyUtils.h"
#include "FSGroupDescriptor.h"

#define DISTANCE_THRESHOLD2   36
#define POINT_RADIUS          4

WidgetGroupPlot::WidgetGroupPlot(QWidget *parent) :
  QWidget(parent), m_bAutoScale(true),
  m_fsgd(NULL),
  m_nCurrentVariableIndex(0),
  m_nPlotType(Point)
{
  setFocusPolicy(Qt::StrongFocus);
  m_ptCurrent = QPointF(-1, -1);
}

WidgetGroupPlot::~WidgetGroupPlot()
{
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
  int nMargin = 10;
  rc_plot.adjust(nMargin, nMargin, -nMargin, -nMargin);
  rc_plot.adjust(15, 15, -25, -27);

  if (m_fsgd->m_classes.isEmpty())
  {
    p.fillRect(rc_plot, Qt::black);
    return;
  }

  if (!m_fsgd->m_title.isEmpty())
  {
    rc_plot.adjust(0, 5, 0, 0);
    p.setPen(Qt::black);
    QRectF rc = rect();
    rc.setBottom(rc_plot.top());
    rc.adjust(0, 5, 0, 0);
    p.drawText(rc, Qt::AlignCenter, m_fsgd->m_title);
  }

  QFont fnt = font();
  fnt.setPixelSize(11);
  int nTextLen = qMax(QFontMetrics(fnt).width(QString::number(m_fsgd->m_variables[0].range[0])),
                      QFontMetrics(fnt).width(QString::number(m_fsgd->m_variables[0].range[1])));
  rc_plot.adjust(nTextLen+6, 0, 0, 0);
  p.fillRect(rc_plot.adjusted(-1, -1, 1, 1), Qt::white);

  p.setPen(Qt::black);
  QMap<QString, FSGDClass> classes;
  for (int i = 0; i < m_fsgd->m_classes.size(); i++)
    classes[m_fsgd->m_classes[i].label] = m_fsgd->m_classes[i];

  int nStart = 0, nEnd = m_fsgd->m_variables.size()-1;
  if (m_nCurrentVariableIndex >= 0)
    nStart = nEnd = m_nCurrentVariableIndex;
  for (int n = nStart; n <= nEnd ; n++)
  {
    double dMin = m_fsgd->m_variables[n].range[0],
           dMax = m_fsgd->m_variables[n].range[1];
    if (m_bAutoScale)
    {
      double val = (dMax-dMin)/5;
      dMax += val;
      dMin -= val;
    }
    if (m_fsgd->m_variables[n].range[0] >= 0 && dMin < 0)
      dMin = 0;
    if (dMin == dMax)
      dMax += 1;
    QList<double> data;
    for (int i = 0; i < m_fsgd->m_data.size(); i++)
      data << m_fsgd->m_data[i].variable_values[n];

    double dSpacing = rc_plot.width() / (data.size()-1+2);
    QPointF* pts = new QPointF[data.size()];
    QString tooltip_strg;
    bool bHit = false;
    for (int i = 0; i < data.size(); i++)
    {
      pts[i] = QPointF(dSpacing + rc_plot.left() + dSpacing*i,
                       rc_plot.bottom() - (data[i]-dMin)/(dMax-dMin)*rc_plot.height());
      QColor c = classes[m_fsgd->m_data[i].class_id].color;
      if (m_nPlotType == Histogram)
      {
        p.setBrush(c);
        QRectF rcf = QRectF(pts[i]-QPointF(dSpacing/2, 0),
                            QPointF(pts[i].x()+dSpacing/2, rc_plot.bottom()-(0-dMin)/(dMax-dMin)*rc_plot.height())).normalized();
        bHit = rcf.contains(m_ptCurrent);
        if (bHit)
          p.setBrush(c.lighter());
        p.drawRect(rcf);
      }
      else if (m_nPlotType == Point)
      {
        if (!bHit)
          bHit = ((pts[i].x()-m_ptCurrent.x())*(pts[i].x()-m_ptCurrent.x()) +
                 (pts[i].y()-m_ptCurrent.y())*(pts[i].y()-m_ptCurrent.y()) < DISTANCE_THRESHOLD2);
        else
          bHit = false;

        DrawMarker(&p, pts[i], classes[m_fsgd->m_data[i].class_id].marker.toLower(),
                   bHit?c.lighter():c);
      }
      // draw current values
      if (bHit)
      {
        tooltip_strg = QString("Subject: %1\nClass: %2")
                       .arg(m_fsgd->m_data[i].subject_id)
                       .arg(m_fsgd->m_data[i].class_id);
        for (int m = n; m <=n; m++)
        {
          tooltip_strg += QString("\n%1: %2").arg(m_fsgd->m_variables[m].label)
                          .arg(m_fsgd->m_data[i].variable_values[m]);
        }
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
        if (qAbs(dMetricPos) < 1e-8)
          strg = "0";
        p.drawText(QRectF(rect().left(), y-10, rc_plot.left()-5-rect().left(), 20),
                   Qt::AlignVCenter | Qt::AlignRight, strg);
      }
      dMetricPos += dMetricStep;
      y = rc_plot.bottom()-(dMetricPos-dMin)/(dMax-dMin)*rc_plot.height();
    }

    p.save();
    p.translate(0, rc_plot.top()+rc_plot.height()/2);
    p.rotate(-90);
    p.drawText(QRect(-100, 0, 200, 20), Qt::AlignCenter, "");
    p.restore();

    // draw X metrics
    if (true)
    {
      nMetricInterval = 50;
      dMetricStep =  (data.size()+1)*m_dTR / (rc_plot.width() / nMetricInterval);
      dMetricStep = MyUtils::RoundToGrid( dMetricStep );
      dMetricPos = 0;
      double x = rc_plot.left();
      while (x < rc_plot.right())
      {
        QString strg = QString::number(dMetricPos);
        p.drawText(QRectF(x-100, rc_plot.bottom()+5, 200, 20),
                   Qt::AlignTop | Qt::AlignHCenter, strg);

        dMetricPos += dMetricStep;
        x = rc_plot.left() + dMetricPos/((data.size()-1)*m_dTR)*rc_plot.width();
      }
      if (dMin < 0 && dMax > 0)
      {
        p.setPen(Qt::black);
        double y = rc_plot.bottom()-(0-dMin)/(dMax-dMin)*rc_plot.height();
        p.drawLine(QPointF(rc_plot.left(), y), QPointF(rc_plot.right(), y));
      }
    }

    p.setPen(Qt::black);
    p.setBrush(Qt::NoBrush);
    p.drawRect(rc_plot);

    QRectF rc = rect().adjusted(0, 0, 0, -3);
    p.drawText(rc, Qt::AlignBottom | Qt::AlignHCenter, "");

    // draw current stats
    if (!tooltip_strg.isEmpty())
    {
      rc = p.boundingRect(QRectF(0,0,1,1), Qt::AlignLeft|Qt::AlignVCenter, tooltip_strg);
      int bw = 6;
      rc.adjust(-bw, -bw, bw, bw);
      rc.moveTo(m_ptCurrent + QPoint(6, 6));
      QPointF pt_c = rc.center();
      qreal w = rc.width()+4, h = rc.height()+4;
      pt_c.rx() = qMin(rc_plot.right()-w/2, qMax(pt_c.rx(), rc_plot.left()+w/2));
      pt_c.ry() = qMin(rc_plot.bottom()-h/2, qMax(pt_c.ry(), rc_plot.top()+h/2));
      rc.moveCenter(pt_c);
      p.setPen(Qt::black);
      p.setBrush(QBrush(QColor(255, 255, 192)));
      p.drawRect(rc);
      rc.adjust(bw, bw, -bw, -bw);
      p.drawText(rc, Qt::AlignLeft | Qt::AlignVCenter, tooltip_strg);
    }
    m_rectPlot = rc_plot;
  }

  // draw class labels
  QPointF pt = rc_plot.topRight()+QPointF(-12, 10);
  QSizeF sz(30, 14);
  if (m_nPlotType == Point)
    sz = QSizeF(14, 14);
  p.setRenderHint(QPainter::Antialiasing, m_nPlotType == Point);
  for (int i = m_fsgd->m_classes.size()-1; i >= 0; i--)
  {
    pt -= QPointF(sz.width(), 0);
    if (m_nPlotType == Histogram)
    {
      p.setPen(Qt::black);
      p.setBrush(m_fsgd->m_classes[i].color);
      p.drawRect(QRectF(pt, sz));
    }
    else if (m_nPlotType == Point)
    {
      DrawMarker(&p, QRectF(pt, sz).center(), m_fsgd->m_classes[i].marker,
                 m_fsgd->m_classes[i].color);
    }
    QRectF rc = p.boundingRect(QRectF(0,0,1,1), Qt::AlignRight|Qt::AlignVCenter, m_fsgd->m_classes[i].label);
    pt -= QPointF(rc.width()+5, 0);
    p.setPen(Qt::black);
    p.drawText(QRectF(pt, rc.size()), Qt::AlignRight|Qt::AlignVCenter, m_fsgd->m_classes[i].label);
    pt -= QPointF(20, 0);
  }
  p.setRenderHint(QPainter::Antialiasing, false);
}

void WidgetGroupPlot::DrawMarker(QPainter *p, const QPointF &pt, const QString &marker, const QColor &c)
{
  p->save();
  double r = POINT_RADIUS;
  p->setPen(QPen(c, 2));
  p->setRenderHint(QPainter::Antialiasing);
  if (marker == "circle")
  {
    p->setBrush(Qt::NoBrush);
    p->drawEllipse(pt, r, r);
  }
  else if (marker == "plus")
  {
    p->drawLine(pt-QPointF(r, 0), pt+QPointF(r, 0));
    p->drawLine(pt-QPointF(0, r), pt+QPointF(0, r));
  }
  else if (marker == "dot" || marker == "point")
  {
    p->setRenderHint(QPainter::Antialiasing);
    p->setBrush(c);
    p->drawEllipse(pt, r-1, r-1);
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
    p->setBrush(Qt::NoBrush);
    p->drawPolygon(pts, 3);
  }
  else if (marker == "rectangle")
  {
    p->setBrush(Qt::NoBrush);
    p->drawRect(QRectF(pt-QPointF(r, r), pt+QPointF(r,r)));
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
    p->setBrush(Qt::NoBrush);
    p->drawPolygon(pts, 4);
  }
  else if (marker == "cross")
  {
    r -= 1;
    p->drawLine(pt-QPointF(r, r), pt+QPointF(r, r));
    p->drawLine(pt-QPointF(r, -r), pt+QPointF(r, -r));
  }
  else if (marker == "asterisk")
  {
    double d = r*0.866;
    p->drawLine(pt-QPointF(d, r/2), pt+QPointF(d, r/2));
    p->drawLine(pt-QPointF(d, -r/2), pt+QPointF(d, -r/2));
    p->drawLine(pt-QPointF(0, r), pt+QPointF(0, r));
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
  if (m_rectPlot.contains(e->posF()))
  {
    m_ptCurrent = e->posF();
    update();
  }
}

void WidgetGroupPlot::mouseMoveEvent(QMouseEvent *e)
{

}

void WidgetGroupPlot::leaveEvent(QEvent *e)
{
  m_ptCurrent = QPointF(-1, -1);
}

void WidgetGroupPlot::keyPressEvent(QKeyEvent *e)
{
}
