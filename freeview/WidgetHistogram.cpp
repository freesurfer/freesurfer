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
#include "WidgetHistogram.h"
#include <QDebug>
#include <QPainter>
#include "MyUtils.h"
#include <QMouseEvent>
#include <QColorDialog>
#include <QMessageBox>
#include <QPalette>

WidgetHistogram::WidgetHistogram(QWidget *parent) :
  QWidget(parent)
{
  m_dInputData = NULL;
  m_nInputSize = 0;
  m_nOutputData = NULL;
  m_bAutoRange = true;
  m_nNumberOfBins = 100;
  m_colorBackground = Qt::white;
  m_colorForeground = Qt::gray;
  m_nMaxCount = 0;
  m_nFixedMaxCount = 0;
  m_nColorTable = NULL;
  m_bMarkerEditable = false;
  m_bUsePercentile = false;
  m_dOutputArea = NULL;

  m_rectGraph = QRect( 50, 20, 100, 100 );
}

WidgetHistogram::~WidgetHistogram()
{
  if ( m_dInputData )
  {
    delete[] m_dInputData;
  }

  if ( m_nOutputData )
  {
    delete[] m_nOutputData;
  }

  if ( m_nColorTable )
  {
    delete[] m_nColorTable;
  }
}


void WidgetHistogram::GetOutputRange( double* dRange )
{
  dRange[0] = m_dOutputRange[0];
  dRange[1] = m_dOutputRange[1];
}

void WidgetHistogram::SetOutputRange( double* dRange )
{
  m_dOutputRange[0] = dRange[0];
  m_dOutputRange[1] = dRange[1];
  m_bAutoRange = false;
  UpdateData();
}

bool WidgetHistogram::GetAutoRange()
{
  return m_bAutoRange;
}

void WidgetHistogram::SetAutoRange( bool bRange )
{
  m_bAutoRange = bRange;
  if ( bRange )
  {
    m_dOutputRange[0] = m_dInputRange[0];
    m_dOutputRange[1] = m_dInputRange[1];
    UpdateData();
  }
}

int WidgetHistogram::GetNumberOfBins()
{
  return m_nNumberOfBins;
}

void WidgetHistogram::SetNumberOfBins( int nBins )
{
  m_nNumberOfBins = nBins;
  UpdateData();
}

void WidgetHistogram::GetOutputData( int* buffer_out )
{
  memcpy( buffer_out, m_nOutputData, m_nNumberOfBins * sizeof( int ) );
}

void WidgetHistogram::UpdateData( bool bRepaint )
{
  if ( m_dInputData == NULL || m_nNumberOfBins < 2 )
  {
    return;
  }

  m_dBinWidth = ( m_dOutputRange[1] - m_dOutputRange[0] ) / m_nNumberOfBins;
  if ( m_nOutputData )
  {
    delete[] m_nOutputData;
  }
  if ( m_dOutputArea)
    delete[] m_dOutputArea;

  m_nOutputData = new int[m_nNumberOfBins];
  m_dOutputArea = new double[m_nNumberOfBins];
  if ( !m_nOutputData || !m_dOutputArea)
  {
    qCritical() << "Can not allocate memory.";
    return;
  }

  // calculate histogram data
  memset( m_nOutputData, 0, m_nNumberOfBins * sizeof( int ) );
  memset( m_dOutputArea, 0, m_nNumberOfBins * sizeof( double ) );
  for ( long i = 0; i < m_nInputSize; i++ )
  {
    int n = (int)( ( m_dInputData[i] - m_dOutputRange[0] ) / m_dBinWidth );
    if ( n >= 0 && n < m_nNumberOfBins )
    {
      m_nOutputData[n] ++;
    }
  }

  m_dOutputTotalArea = 0;
  for (int i = 0; i < m_nNumberOfBins; i++)
  {
    m_dOutputTotalArea += m_nOutputData[i];
    m_dOutputArea[i] = m_dOutputTotalArea;
  }

  // find max and second max
  m_nMaxCount = 0;
  int nSecondMax = 0;
  for ( int i = 0; i < m_nNumberOfBins; i++ )
  {
    if ( m_nMaxCount < m_nOutputData[i] )
    {
      m_nMaxCount  = m_nOutputData[i];
    }
    else if ( nSecondMax < m_nOutputData[i] )
    {
      nSecondMax = m_nOutputData[i];
    }
  }

  if ( m_nMaxCount > nSecondMax * 5 )
  {
    m_nMaxCount = nSecondMax;
  }

  // allocate color table
  if ( m_nColorTable )
  {
    delete[] m_nColorTable;
  }

  m_nColorTable = new unsigned char[ m_nNumberOfBins*4 ];
  UpdateColorTable();

  if ( bRepaint )
  {
    this->repaint();
  }
}

void WidgetHistogram::UpdateColorTable()
{
  for ( int i = 0; i < m_nNumberOfBins; i++ )
  {
    m_nColorTable[i*4] = m_colorForeground.red();
    m_nColorTable[i*4+1] = m_colorForeground.green();
    m_nColorTable[i*4+2] = m_colorForeground.blue();
    m_nColorTable[i*4+3] = 255;
  }
}

void WidgetHistogram::paintEvent(QPaintEvent* event)
{
  Q_UNUSED(event);
  //    QWidget::paintEvent( event );

  QPainter painter( this );
  QRect rc = rect();
  if ( m_nOutputData )
  {
    int nMaxCnt = (m_nFixedMaxCount > 0 ? m_nFixedMaxCount : m_nMaxCount);
    int x, y;
    int nOrigin[2] = { m_rectGraph.x(), m_rectGraph.y() };
    int nCavWidth = rc.width() - m_rectGraph.x() - 20;
    int nCavHeight = rc.height() - m_rectGraph.y() - 20;
    m_rectGraph.setSize( QSize(nCavWidth, nCavHeight) );

    // draw background
    painter.setBrush( QBrush( m_colorBackground ) );
    painter.setPen( QPen(m_colorBackground) );
    painter.drawRect( m_rectGraph );

    // draw y metrics
    int nMetricInterval = 25;
    QPalette pal = palette();
    double dMetricStep = ((double)nMaxCnt) / ( nCavHeight / nMetricInterval );
    dMetricStep = MyUtils::RoundToGrid( dMetricStep );
    double dMetricStart = 0;
    y = m_rectGraph.bottom();
    painter.setPen( QPen(pal.color(QPalette::WindowText)) );
    while ( y > m_rectGraph.top() && dMetricStep > 0 )
    {
      if( y < m_rectGraph.bottom() )
      {
        painter.drawLine( m_rectGraph.left(), y, m_rectGraph.left()-4, y );

        QPen oldPen = painter.pen();
        QPen newpen( Qt::gray );
        newpen.setStyle(Qt::DotLine);
        painter.setPen( newpen );
        painter.drawLine( m_rectGraph.left(), y, m_rectGraph.right(), y );
        painter.setPen( oldPen );
      }
      QString value_strg = QString::number( dMetricStart );
      QRect tmp_rc = painter.boundingRect(QRect(), Qt::AlignRight, value_strg );
      tmp_rc.moveCenter( QPoint(m_rectGraph.left()-tmp_rc.width()/2-7, y) );
      painter.drawText( tmp_rc, value_strg );

      dMetricStart += dMetricStep;
      y =  (int)(m_rectGraph.bottom() - dMetricStart / nMaxCnt *nCavHeight );
    }

    // draw bars
    double dStepWidth = ( (double) nCavWidth) / m_nNumberOfBins;
    x = nOrigin[0];
    int nLastPos = x;
    for ( int i = 0; i < m_nNumberOfBins; i++ )
    {
      painter.setPen( QPen( QColor( m_nColorTable[i*4],  m_nColorTable[i*4+1], m_nColorTable[i*4+2] ) ) );
      painter.setBrush(QBrush( QColor( m_nColorTable[i*4],  m_nColorTable[i*4+1], m_nColorTable[i*4+2] ) ) );
      y = (int)( nOrigin[1] + nCavHeight * ( 1.0 - (double)m_nOutputData[i] / nMaxCnt ) );
      int h = (int)( (double)m_nOutputData[i] / nMaxCnt * nCavHeight );
      if ( y < nOrigin[1] )
      {
        y = nOrigin[1];
        h = nCavHeight;
      }
      int x2 = (int)( nOrigin[0] + dStepWidth * (i+1) );

      painter.drawRect( x, y, x2-x, h );

      if (m_bUsePercentile)
      {
        bool bDraw = false;
        double dVal = 0;
        int nPos = x;
        if (i == 0 || i == m_nNumberOfBins-1)
        {
            bDraw = true;
            if (i == m_nNumberOfBins-1)
            {
              dVal = 100;
              nPos = m_rectGraph.right();
            }
        }
        else
        {
            double val1 = 100.0*m_dOutputArea[i-1]/m_dOutputTotalArea,
                val2 = 100.0*m_dOutputArea[i]/m_dOutputTotalArea;
            if (((int)val1) != ((int)val2) && x-nLastPos > 40)
            {
              nLastPos = x;
              bDraw = true;
              dVal = (int)val2;
              nPos = (x2-x)*(((int)val2)-val1)/(val2-val1) + x - (x2-x)/2.0;
            }
        }
        if (bDraw)
        {
          painter.setPen(pal.color(QPalette::WindowText));
          if (i > 0 && i < m_nNumberOfBins-1)
            painter.drawLine( nPos, m_rectGraph.bottom(), nPos, m_rectGraph.bottom()+4 );
          QString value_strg = QString::number( (int)dVal );
          QRect tmp_rc = painter.boundingRect( QRect(), Qt::AlignCenter, value_strg );
          tmp_rc.moveCenter( QPoint((x+x2)/2, m_rectGraph.bottom()+5+tmp_rc.height()/2));
          painter.drawText( tmp_rc, value_strg );
        }
      }
      x = x2;
    }

    // draw axis
    painter.setBrush( QBrush(Qt::NoBrush) );
    painter.setPen( QPen(pal.color(QPalette::WindowText)) );
    painter.drawRect( m_rectGraph );

    // draw zero line
    LineMarker lm;
    lm.position = 0;
    lm.movable = false;
    lm.style = Qt::DashLine;
    lm.color = Qt::black;
    DrawMarker(&painter, lm);

    // draw markers
    for ( int i = 0; i < m_markers.size(); i++ )
    {
      LineMarker lm = m_markers[i];
      if (m_bSymmetricMarkers && lm.position < 0)
      {
        lm.position = -lm.position;
        m_markers[i] = lm;
      }
      DrawMarker(&painter, lm);
      if ( m_bSymmetricMarkers )
      {
        lm.position = -lm.position;
        DrawMarker(&painter, lm);
      }
    }

    // draw x metrics
    if (!m_bUsePercentile)
    {
      painter.setPen(QPen(pal.color(QPalette::WindowText)));
      nMetricInterval = 50;
      dMetricStep = ( m_dOutputRange[1] - m_dOutputRange[0] ) / ( nCavWidth / nMetricInterval );
      dMetricStep = MyUtils::RoundToGrid( dMetricStep );
      dMetricStart = (int)( ( m_dOutputRange[0] / dMetricStep ) ) * dMetricStep;
      if ( m_dOutputRange[0] < 0 )
      {
        dMetricStart -= dMetricStep;
      }
      x = ( int )( nOrigin[0] + ( dMetricStart - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] )
          *nCavWidth );
      while ( x <= m_rectGraph.right() && dMetricStep > 0 )
      {
        if( x >= m_rectGraph.left() )
        {
          painter.drawLine( x, m_rectGraph.bottom(), x, m_rectGraph.bottom()+4 );
        }
        QString value_strg = QString::number( dMetricStart );
        QRect tmp_rc = painter.boundingRect( QRect(), Qt::AlignCenter, value_strg );
        if ( x - tmp_rc.width() / 2 > 0 )
        {
          tmp_rc.moveCenter( QPoint(x, m_rectGraph.bottom()+5+tmp_rc.height()/2));
          painter.drawText( tmp_rc, value_strg );
        }

        dMetricStart += dMetricStep;
        if ( fabs( dMetricStart ) < 1e-10 )
        {
          dMetricStart = 0;
        }

        x = ( int )(m_rectGraph.left() + ( dMetricStart - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] )*nCavWidth );
      }
    }
  }
}

void WidgetHistogram::DrawMarker(QPainter *p, LineMarker marker)
{
  if (marker.movable)
  {
    p->setPen( QPen( QBrush(Qt::darkGray), 1, marker.style ) );
  }
  else
  {
    p->setPen( QPen( QBrush(marker.color), 1, marker.style) );
  }
  int x = ( int )( m_rectGraph.left() +
                   m_rectGraph.width() * ( marker.position - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] ) );
  int y = m_rectGraph.y();
  if ( x >= m_rectGraph.left() && x <= m_rectGraph.right()+1 )
  {
    p->drawLine( x, y, x, y + m_rectGraph.height() );
    if ( marker.movable)
    {
      DrawMarkerThumb(p, marker, x, y);
    }
  }
}

void WidgetHistogram::DrawMarkerThumb(QPainter* p, LineMarker marker, int x, int y)
{
  p->setPen(QPen(Qt::black));
  p->setBrush(QBrush(marker.color));
  p->drawPolygon(MakeMarkerThumb(x, y));
}

void WidgetHistogram::SetForegroundColor( const QColor& color )
{
  m_colorForeground = color;
}

void WidgetHistogram::SetColorTableData( unsigned char* colortable, bool bRefresh )
{
  memcpy( m_nColorTable, colortable, m_nNumberOfBins * 4 );
  if ( bRefresh )
  {
    repaint();
  }
}

void WidgetHistogram::SetMarkers( const LineMarkers& markers, bool bRefresh  )
{
  m_markers = markers;

  if ( bRefresh )
  {
    repaint();
  }
}

bool WidgetHistogram::FindMarker(int px, int py, int* nIndex, bool* bMirrored)
{
  for ( int i = 0; i < m_markers.size(); i++ )
  {
    LineMarker marker = m_markers[i];
    if (marker.movable)
    {
      if (MarkerHit(px, py, marker))
      {
        *nIndex = i;
        if (bMirrored)
        {
          *bMirrored = false;
        }
        return true;
      }
      else
      {
        marker.position = -marker.position;
        if (MarkerHit(px, py, marker))
        {
          *nIndex = i;
          if (bMirrored)
          {
            *bMirrored = true;
          }
          return true;
        }
      }
    }
  }
  return false;
}

bool WidgetHistogram::MarkerHit(int px, int py, LineMarker marker)
{
  int x = ( int )( m_rectGraph.left() +
                   m_rectGraph.width() * ( marker.position - m_dOutputRange[0] ) / ( m_dOutputRange[1] - m_dOutputRange[0] ) );
  int y = m_rectGraph.y();
  QPolygon polygon = MakeMarkerThumb(x, y);
  return polygon.containsPoint(QPoint(px, py), Qt::OddEvenFill);
}

QPolygon WidgetHistogram::MakeMarkerThumb(int x, int y)
{
  QVector<QPoint> points;
  int r = 4;
  points << QPoint(x-r, y-r)
         << QPoint(x+r, y-r)
         << QPoint(x+r, y)
         << QPoint(x, y+r)
         << QPoint(x-r, y);
  return QPolygon(points);
}

void WidgetHistogram::mousePressEvent(QMouseEvent* event )
{
  /*
  if ( m_rectGraph.contains( event->x(), event->y() ) )
  {
  double dValue = ( event->x() - m_rectGraph.left() ) *
                  ( m_dOutputRange[1] - m_dOutputRange[0] ) / m_rectGraph.width() + m_dOutputRange[0];
  if ( event->button() == Qt::LeftButton )
  {
      emit MouseButtonPressed(Qt::LeftButton, dValue);
  }
  else if ( event->button() == Qt::MiddleButton )
  {
      emit MouseButtonPressed(Qt::MiddleButton, dValue);
  }
  else if ( event->button() == Qt::RightButton )
  {
      emit MouseButtonPressed(Qt::RightButton, dValue);
  }
  }
  */
  if (event->button() == Qt::LeftButton)
  {
    if (FindMarker(event->x(), event->y(), &m_nActiveMarker, &m_bActiveMarkerMirrored))
    {
      if (m_bMarkerEditable && (event->modifiers() & Qt::ShiftModifier) )
      {
        if ( m_markers.size() == 1)
        {
          QMessageBox::warning(this, "Error", "Can not remove the only color stop left.");
        }
        else
        {
          m_markers.remove(m_nActiveMarker);
          m_nActiveMarker = -1;
          emit MarkerChanged();
        }
      }
    }
    else
    {
      m_nActiveMarker = -1;
      if (m_markers.size() > 0)
      {
        int nMarker = 0;
        int x = event->x();
        double pos = (x - m_rectGraph.left())*(m_dOutputRange[1] - m_dOutputRange[0])/m_rectGraph.width() + m_dOutputRange[0];
        if (pos < m_dOutputRange[0])
        {
          pos = m_dOutputRange[0];
        }
        else if (pos > m_dOutputRange[1])
        {
          pos = m_dOutputRange[1];
        }
        if (m_bActiveMarkerMirrored)
        {
          m_markers[nMarker].position = -pos;
        }
        else
        {
          m_markers[nMarker].position = pos;
        }

        this->repaint();
        emit MarkerChanged();
      }
    }
  }
  else if (event->button() == Qt::RightButton)
  {
    if (m_markers.size() > 1)
    {
      int nMarker = m_markers.size()-1;
      int x = event->x();
      double pos = (x - m_rectGraph.left())*(m_dOutputRange[1] - m_dOutputRange[0])/m_rectGraph.width() + m_dOutputRange[0];
      if (pos < m_dOutputRange[0])
      {
        pos = m_dOutputRange[0];
      }
      else if (pos > m_dOutputRange[1])
      {
        pos = m_dOutputRange[1];
      }
      if (m_bActiveMarkerMirrored)
      {
        m_markers[nMarker].position = -pos;
      }
      else
      {
        m_markers[nMarker].position = pos;
      }

      this->repaint();
      emit MarkerChanged();
    }
  }
}


void WidgetHistogram::mouseDoubleClickEvent(QMouseEvent * event)
{
  int n;
  if (m_bMarkerEditable && event->button() == Qt::LeftButton &&
      FindMarker(event->x(), event->y(), &n) )
  {
    QColor c = QColorDialog::getColor(m_markers[n].color);
    if (c.isValid())
    {
      m_markers[n].color = c;
      this->update();
      emit MarkerChanged();
    }
  }
}

void WidgetHistogram::mouseMoveEvent(QMouseEvent* event)
{
  if (m_nActiveMarker >= 0)
  {
    int x = event->x();
    double pos = (x - m_rectGraph.left())*(m_dOutputRange[1] - m_dOutputRange[0])/m_rectGraph.width() + m_dOutputRange[0];
    if (pos < m_dOutputRange[0])
    {
      pos = m_dOutputRange[0];
    }
    else if (pos > m_dOutputRange[1])
    {
      pos = m_dOutputRange[1];
    }
    if (m_bActiveMarkerMirrored)
    {
      m_markers[m_nActiveMarker].position = -pos;
    }
    else
    {
      m_markers[m_nActiveMarker].position = pos;
    }

    this->repaint();
    emit MarkerChanged();
  }
}

void WidgetHistogram::mouseReleaseEvent(QMouseEvent * event)
{
  Q_UNUSED(event);
  m_nActiveMarker = -1;
}

void WidgetHistogram::AddMarker(double pos, const QColor &color)
{
  bool inserted = false;
  for (int i = 0; i < m_markers.size(); i++ )
  {
    if (m_markers[i].position > pos )
    {
      LineMarker m;
      m.color = color;
      m.position = pos;
      m_markers.insert(i, m);
      inserted = true;
      break;
    }
  }
  if ( !inserted )
  {
    LineMarker m;
    m.color = color;
    m.position = pos;
    m_markers << m;
  }
  this->repaint();
  emit MarkerChanged();
}

void WidgetHistogram::FlipMarkers()
{
  LineMarkers markers;
  double dMin = m_markers[0].position;
  double dMax = dMin;
  for (int i = 0; i < m_markers.size(); i++)
  {
    if ( dMin > m_markers[i].position)
    {
      dMin = m_markers[i].position;
    }
    else if (dMax < m_markers[i].position)
    {
      dMax = m_markers[i].position;
    }
  }
  for (int i = 0; i < m_markers.size(); i++)
  {
    LineMarker m = m_markers[i];
    m.position = dMax - (m.position-dMin);
    markers << m;
  }
  m_markers = markers;
  this->repaint();
  emit MarkerChanged();
}

void WidgetHistogram::SetMarkerEditable(bool bFlag)
{
  m_bMarkerEditable = bFlag;
  setToolTip(bFlag?"Double-click on the sliding markers to edit.\r\nShift+Click to delete":"");
}

void WidgetHistogram::SetFixedMaxCount(int cnt)
{
  m_nFixedMaxCount = cnt;
  this->repaint();
}
