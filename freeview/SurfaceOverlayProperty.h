/**
 * @file  SurfaceOverlayProperty.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/04/29 17:27:01 $
 *    $Revision: 1.5 $
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
 *
 */

#ifndef SurfaceOverlayProperty_h
#define SurfaceOverlayProperty_h

#include "vtkSmartPointer.h"
#include <QObject>
#include <QColor>
#include <QPair>
#include <QVector>

#ifndef QGradientStop
typedef QPair<qreal, QColor> QGradientStop;
typedef QVector<QGradientStop> QGradientStops;
#endif

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class SurfaceOverlay;

class SurfaceOverlayProperty  : public QObject
{
  Q_OBJECT
public:
  SurfaceOverlayProperty ( SurfaceOverlay* overlay );
  ~SurfaceOverlayProperty ();

  enum COLOR_SCALE  { CS_Heat = 0, CS_GreenRed, CS_BlueRed, CS_ColorWheel, CS_Custom };
  enum COLOR_METHOD { CM_Linear = 0, CM_LinearOpaque, CM_Piecewise };

  double GetOpacity() const;
  void SetOpacity( double opacity );

  int GetColorScale() const;
  void SetColorScale( int nScale );

  void SetSurfaceOverlay( SurfaceOverlay* overlay );

  void SetMinPoint( double dValue );
  double GetMinPoint();

  void SetMidPoint( double dValue );
  double GetMidPoint();

  void SetMaxPoint( double dValue );
  double GetMaxPoint();

  int GetColorMethod();
  void SetColorMethod( int n );

  bool GetColorInverse();
  void SetColorInverse( bool bInverse );

  bool GetColorTruncate();
  void SetColorTruncate( bool bTruncate );

  bool GetClearLower()
  {
    return m_bClearLower;
  }
  void SetClearLower( bool bClear );

  bool GetClearHigher()
  {
    return m_bClearHigher;
  }
  void SetClearHigher(bool bClear);

  QGradientStops GetCustomColorScale()
  {
    return m_customScale;
  }
  void SetCustomColorScale(QGradientStops stops);

  bool GetSmooth()
  {
    return m_bSmooth;
  }
  void SetSmooth(bool bSmooth);

  int GetSmoothSteps()
  {
    return m_nSmoothSteps;
  }
  void SetSmoothSteps(int n);

//  void MapOverlayColor( unsigned char* colordata, int nPoints );
  void MapOverlayColor( float* data, unsigned char* colordata, int nPoints );
  void MapOverlayColorSymmetric( float* data, unsigned char* colordata, int nPoints );
  void MapOverlayColorFullScale( float* data, unsigned char* colordata, int nPoints );

  vtkRGBAColorTransferFunction* GetLookupTable()
  {
    return m_lut;
  }

  void Reset();

  void EmitColorMapChanged()
  {
    emit ColorMapChanged();
  }

Q_SIGNALS:
  void ColorMapChanged();
  void SmoothChanged();

private:

  double      m_dOpacity;
  int         m_nColorScale;
  int         m_nColorMethod;
  double      m_dMinPoint;
  double      m_dMidPoint;
  double      m_dMaxPoint;
  int         m_colorMin[3];
  int         m_colorMid[3];
  int         m_colorMax[3];
  bool        m_bColorInverse;
  bool        m_bColorTruncate;
  QGradientStops m_customScale;
  double      m_dMinStop;
  double      m_dMaxStop;
  bool        m_bClearLower;
  bool        m_bClearHigher;
  bool        m_bSmooth;
  int         m_nSmoothSteps;

  vtkRGBAColorTransferFunction* m_lut;

  SurfaceOverlay* m_overlay;
};

#endif
