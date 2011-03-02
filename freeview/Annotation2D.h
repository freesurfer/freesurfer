/**
 * @file  Annotation2D.h
 * @brief Annotation for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.10 $
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

#ifndef Annotation2D_h
#define Annotation2D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"

class vtkTextActor;
class vtkRenderer;
class vtkAxisActor2D;
class vtkPropCollection;

class Annotation2D
{
public:
  Annotation2D();
  virtual ~Annotation2D();

  void Update( vtkRenderer* renderer, int nPlane );

  void AppendAnnotations( vtkRenderer* renderer );

  void ShowScaleLine( bool bShow );

  bool GetShowScaleLine();

  void Show( bool bShow );
  bool IsVisible();

private:
  void UpdateScaleActors( double length, int nNumOfTicks, const char* title );

  vtkSmartPointer<vtkTextActor> m_actorCoordinates[6];
  vtkSmartPointer<vtkActor2D>  m_actorScaleLine;
  vtkSmartPointer<vtkTextActor> m_actorScaleTitle;
  vtkSmartPointer<vtkPropCollection> m_actorsAll;
};

#endif


