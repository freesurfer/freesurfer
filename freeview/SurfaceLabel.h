/**
 * @file  SurfaceLabel.h
 * @brief The common properties available to surface label
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:18 $
 *    $Revision: 1.4 $
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

#ifndef SurfaceLabel_h
#define SurfaceLabel_h

#include <QObject>

extern "C"
{
#include "label.h"
}

class LayerSurface;

class SurfaceLabel  : public QObject
{
  Q_OBJECT
public:
  SurfaceLabel ( LayerSurface* surf );
  ~SurfaceLabel ();

//  void SetSurface( LayerSurface* surf );

  QString GetName();

  void SetName( const QString& name );

  bool LoadLabel( const QString& filename );

  void SetColor( double r, double g, double b );

  double* GetColor()
  {
    return m_rgbColor;
  }

  void MapLabel( unsigned char* colordata, int nVertexCount );

Q_SIGNALS:
  void SurfaceLabelChanged();

private:
  LABEL*        m_label;
  QString       m_strName;
  LayerSurface* m_surface;
  double        m_rgbColor[3];
  bool          m_bTkReg;
};

#endif
