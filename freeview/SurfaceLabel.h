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
 *    $Author: rpwang $
 *    $Date: 2013/11/05 20:25:29 $
 *    $Revision: 1.8 $
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
#include <vtkSmartPointer.h>

extern "C"
{
#include "label.h"
}

class LayerSurface;
class vtkActor;
class vtkPolyData;

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

  bool GetShowOutline()
  {
    return m_bShowOutline;
  }

  void SetShowOutline(bool bOutline);

  vtkActor* GetOutlineActor();

  bool IsVisible()
  {
    return m_bVisible;
  }

Q_SIGNALS:
  void SurfaceLabelChanged();

public slots:
  void SetVisible(bool flag);

private:
  QList<int> DoConnectEdgeVertices(const QList<int>& indices_in, const QList<int>& vertices);
  vtkPolyData* MakeEdgePolyData(const QList<int>& indices_in, const QList<int>& vertices);

  LABEL*        m_label;
  QString       m_strName;
  LayerSurface* m_surface;
  double        m_rgbColor[3];
  bool          m_bTkReg;
  bool          m_bShowOutline;
  bool          m_bVisible;
  vtkSmartPointer<vtkActor> m_actorOutline;
};

#endif
