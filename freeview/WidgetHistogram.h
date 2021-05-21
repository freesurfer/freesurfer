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
#ifndef WIDGETHISTOGRAM_H
#define WIDGETHISTOGRAM_H

#include <QWidget>
#include <QVector>
#include <QColor>
#include <QPolygon>

struct LineMarker;
typedef QVector<LineMarker> LineMarkers;

struct LineMarker
{
  LineMarker() : style(Qt::DashLine), movable(true)
  {
  }

  double        position;
  QColor        color;
  Qt::PenStyle  style;
  bool          movable;
};

class WidgetHistogram : public QWidget
{
  Q_OBJECT
public:
  explicit WidgetHistogram(QWidget *parent = 0);
  ~WidgetHistogram();

  template <class T> void SetInputData( T* data, long size, double* range = NULL);

  void GetOutputRange( double* dRange );

  void SetOutputRange( double* dRange );

  bool GetAutoRange();


  int GetNumberOfBins();

  void SetNumberOfBins( int nBins );

  int GetOutputSize();

  void GetOutputData( int* buffer_out );

  int GetMaximumCount()
  {
    return m_nMaxCount;
  }


  void SetMarkers( const LineMarkers& markers, bool bRefresh = true );

  void SetSymmetricMarkers(bool bFlag)
  {
    m_bSymmetricMarkers = bFlag;
  }

  void SetMarkerEditable(bool bFlag);

  void SetColorTableData( unsigned char* colortable, bool bRefresh = true );

  void AddMarker(double pos, const QColor& color);

  void SetFixedMaxCount(int cnt);

  LineMarkers GetMarkers()
  {
    return m_markers;
  }

  void GetInputRange(double* range)
  {
    range[0] = m_dInputRange[0];
    range[1] = m_dInputRange[1];
  }

signals:
  void MouseButtonPressed(int button, double value);
  void MarkerChanged();

public slots:
  void SetAutoRange( bool bRange );
  void SetForegroundColor( const QColor& color );
  void FlipMarkers();
  void SetUsePercentile(bool b)
  {
    m_bUsePercentile = b;
    update();
  }

protected:
  void Initialize();
  virtual void paintEvent(QPaintEvent *);
  virtual void mousePressEvent(QMouseEvent *);
  virtual void mouseMoveEvent(QMouseEvent *);
  virtual void mouseReleaseEvent(QMouseEvent *);
  virtual void mouseDoubleClickEvent(QMouseEvent *);

  void UpdateData( bool bRepaint = true );
  void UpdateColorTable( );
  void DrawMarker(QPainter* p, LineMarker marker);
  void DrawMarkerThumb(QPainter* p, LineMarker marker, int x, int y);
  bool FindMarker(int x, int y, int* nIndex, bool* bMirrored = 0);
  bool MarkerHit(int x, int y, LineMarker marker);
  QPolygon MakeMarkerThumb(int x, int y);

  double*     m_dInputData;
  long        m_nInputSize;
  double      m_dInputRange[2];

  int*        m_nOutputData;
  double      m_dOutputRange[2];
  bool        m_bAutoRange;
  int         m_nNumberOfBins;
  double      m_dOutputTotalArea;
  double*     m_dOutputArea;
  unsigned char* m_nColorTable;       // color table for histogram drawing as RGBA

  double      m_dBinWidth;
  int         m_nMaxCount;
  int         m_nFixedMaxCount;

  QColor    m_colorBackground;
  QColor    m_colorForeground;

  QRect      m_rectGraph;

  LineMarkers m_markers;
  bool        m_bSymmetricMarkers;
  bool        m_bMarkerEditable;
  int         m_nActiveMarker;
  bool        m_bActiveMarkerMirrored;
  bool        m_bUsePercentile;
};

template <class T> void WidgetHistogram::SetInputData( T* data, long size, double* range)
{
  if ( m_dInputData )
  {
    delete[] m_dInputData;
  }

  m_dInputData = new double[size];
  if ( m_dInputData )
  {
    m_dInputRange[0] = m_dInputRange[1] = data[0];
    for ( long i = 0; i < size; i++ )
    {
      m_dInputData[i] = data[i];
      if ( m_dInputRange[0] > m_dInputData[i] )
      {
        m_dInputRange[0] = m_dInputData[i];
      }
      else if ( m_dInputRange[1] < m_dInputData[i] )
      {
        m_dInputRange[1] = m_dInputData[i];
      }
    }

    m_nInputSize = size;
    m_dOutputRange[0] = m_dInputRange[0];
    m_dOutputRange[1] = m_dInputRange[1];
  }

  if (range)
  {
    m_dOutputRange[0] = range[0];
    m_dOutputRange[1] = range[1];
  }

  UpdateData();
}

#endif // WIDGETHISTOGRAM_H
