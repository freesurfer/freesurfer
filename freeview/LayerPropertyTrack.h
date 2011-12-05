/**
 * @file  LayerPropertyTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/05 20:03:33 $
 *    $Revision: 1.6 $
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
#ifndef LAYERPROPERTYTRACK_H
#define LAYERPROPERTYTRACK_H

#include "LayerProperty.h"
#include <QColor>

class LayerPropertyTrack : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyTrack(QObject* parent = 0);

  enum ColorCode { Directional = 0, SolidColor };
  enum DirectionScheme { EndPoints = 0, MidSegment, EverySegment };
  enum DirectionMapping  { RGB = 0, RBG, GRB, GBR, BRG, BGR };
  enum RenderRep { Line = 0, Tube };

  int GetColorCode()
  {
    return m_nColorCode;
  }

  int GetDirectionScheme()
  {
    return m_nDirectionScheme;
  }

  int GetDirectionMapping()
  {
    return m_nDirectionMapping;
  }

  QColor GetSolidColor()
  {
    return m_color;
  }

  int GetRenderRep()
  {
    return m_nRenderRep;
  }

  double GetTubeRadius()
  {
    return m_dTubeRadius;
  }

  int GetNumberOfSides()
  {
    return m_nNumberOfSides;
  }

signals:
  void ColorCodeChanged(int);
  void DirectionSchemeChanged(int);
  void DirectionMappingChanged(int);
  void SolidColorChanged(const QColor& c);
  void RenderRepChanged();

public slots:
  void SetColorCode(int nCode);
  void SetDirectionScheme(int nVal);
  void SetDirectionMapping(int nVal);
  void SetSolidColor(const QColor& c);
  void SetRenderRep(int nVal);
  void SetTubeRadius(double dVal);
  void SetNumberOfSides(int nVal);

private:
  int     m_nColorCode;
  int     m_nDirectionScheme;
  int     m_nDirectionMapping;
  int     m_nRenderRep;
  double  m_dTubeRadius;
  int     m_nNumberOfSides;
  QColor  m_color;
};

#endif // LAYERPROPERTYTRACK_H
