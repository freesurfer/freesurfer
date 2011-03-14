/**
 * @file  Cursor3D.h
 * @brief Cursor for 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:57 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2008-2009,
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

#ifndef Cursor3D_h
#define Cursor3D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"
#include <QColor>
#include <QObject>

class vtkRenderer;
class vtkActor;
class RenderView3D;

class Cursor3D : public QObject
{
    Q_OBJECT
public:
  Cursor3D( RenderView3D* view );
  virtual ~Cursor3D();

  void SetPosition( double* pos );

  double* GetPosition();
  void GetPosition( double* pos );

  void GetColor( double* rgb );
  void SetColor( double r, double g, double b );

  QColor GetColor();

  void Update();

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );

  bool IsShown();

public slots:
  void SetColor( const QColor& color );

Q_SIGNALS:
  void Updated();

private:
  void RebuildActor();

  vtkSmartPointer<vtkActor> m_actorCursor;

  RenderView3D* m_view;

  double  m_dPosition[3];
};

#endif


