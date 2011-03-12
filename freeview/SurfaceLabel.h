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
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:53 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
