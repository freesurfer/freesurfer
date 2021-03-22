/**
 * @brief Annotation for 2D view.
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
 *
 */

#ifndef Annotation2D_h
#define Annotation2D_h

#include <QObject>
#include <QColor>
#include "vtkSmartPointer.h"

class vtkTextActor;
class vtkActor2D;
class vtkRenderer;
class vtkAxisActor2D;
class vtkPropCollection;

class Annotation2D : public QObject
{
  Q_OBJECT
public:
  Annotation2D( QObject* parent );
  virtual ~Annotation2D();

  void Update( vtkRenderer* renderer, int nPlane );

  void AppendAnnotations( vtkRenderer* renderer, bool bScaleBar = false );

  bool GetShowScaleLine();

  void Show( bool bShow );
  bool IsVisible();

  QColor GetColor();

  void SetAutoScaleText(bool b);

  void SetTextSize(int nsize);

public slots:
  void ShowScaleLine( bool bShow );
  void SetColor(const QColor& c);

Q_SIGNALS:
  void Updated();

private:
  void UpdateScaleActors( double length, int nNumOfTicks, const char* title );

  vtkSmartPointer<vtkTextActor> m_actorCoordinates[6];
  vtkSmartPointer<vtkActor2D>  m_actorScaleLine;
  vtkSmartPointer<vtkTextActor> m_actorScaleTitle;
  vtkSmartPointer<vtkPropCollection> m_actorsAll;
};

#endif


