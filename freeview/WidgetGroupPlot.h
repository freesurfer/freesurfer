/**
 * @brief Widget drawing group plot
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

#ifndef WidgetGroupPlot_H
#define WidgetGroupPlot_H

#include <QWidget>
#include <QList>
#include <QColor>

class FSGroupDescriptor;

class WidgetGroupPlot : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetGroupPlot(QWidget* parent = 0);
  ~WidgetGroupPlot();

  void paintEvent(QPaintEvent * e);

  void SetFsgdData(FSGroupDescriptor* fsgd);

  void mousePressEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);
  void leaveEvent(QEvent *e);

  enum PlotType { Point = 0, Histogram };

public slots:
  void SetAutoScale(bool bAuto);
  void SetCurrentVariableIndex(int n);
  void SetPlotType(int n);
  void SetCurrentVertex(int nVertex);
  void SetCurrentDataIndex(int n)
  {
      m_nCurrentDataIndex = n;
      update();
  }

signals:
  void FrameChanged(int frame);
  void CurrentDataIndexChanged(int nIndex);

public:
  static void DrawMarker(QPainter* p, const QPointF& pt, const QString& marker, const QColor& c,
                         double r, bool bHighlight = false);
  void DrawMarker(QPainter* p, const QPointF& pt, const QString& marker, const QColor& c, bool bHighlight);

  double          m_dTR;
  double          m_dMin;
  double          m_dMax;
  bool            m_bAutoScale;
  QRectF          m_rectPlot;
  FSGroupDescriptor*  m_fsgd;
  int       m_nCurrentDataIndex;

  int       m_nCurrentVariableIndex;
  int       m_nPlotType;
};

#endif // WidgetGroupPlot_H
