/**
 * @brief Region2DLine data object.
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

#ifndef Region2DLine_h
#define Region2DLine_h

#include "Region2D.h"
#include "vtkSmartPointer.h"

class vtkActor2D;
class RenderView2D;
class vtkRenderer;
class vtkTextActor;

class Region2DLine : public Region2D
{
public:
  Region2DLine( RenderView2D* view );
  virtual ~Region2DLine();

  void SetLine( int x1, int y1, int x2, int y2 );

  void SetPoint1( int x1, int y1 );
  void SetPoint2( int x2, int y2 );

  void Offset( int x, int y );

  bool Contains( int x, int y, int* nIndexOut = NULL );
  void UpdatePoint( int nIndex, int nX, int nY );

  void AppendProp( vtkRenderer* renderer );

  void Show( bool bshow = true );
  void Highlight( bool bHighlight = true );

  void Update();
  void UpdateStats();

  void UpdateSlicePosition( int nPlane, double pos );

  void GetWorldPoint( int nIndex, double* pt );

  QString DataToString();

  Region2D* ObjectFromString(RenderView2D* view, const QString& text);

protected:
  void UpdateWorldCoords();

  vtkSmartPointer<vtkActor2D>   m_actorLine;
  int       m_nX1, m_nX2, m_nY1, m_nY2;
  double    m_dPt1[3];
  double    m_dPt2[3];

};

#endif


