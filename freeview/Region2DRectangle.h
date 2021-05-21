/**
 * @brief Region2DRectangle data object.
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

#ifndef Region2DRectangle_h
#define Region2DRectangle_h

#include "Region2D.h"
#include <vector>
#include "vtkSmartPointer.h"

class vtkTextActor;
class vtkActor2D;
class RenderView2D;
class vtkRenderer;

class Region2DRectangle : public Region2D
{
public:
  Region2DRectangle( RenderView2D* view );
  virtual ~Region2DRectangle();

  void UpdatePoint( int nIndex, int nX, int nY );

  void Offset( int nX, int nY );

  bool Contains( int nX, int nY, int* nIndexOut = NULL );

  void SetRect( int left, int top, int w, int h );

  void SetTopLeft( int left, int top );
  void SetBottomRight( int right, int bottom );

  void AppendProp( vtkRenderer* renderer );

  void Show( bool bshow );

  void Highlight( bool bHighlight );

  void Update();

  void UpdateStats();

  void UpdateSlicePosition( int nPlane, double pos );

  void GetWorldPoint( int nIndex, double* pt );

  void SetEnableStats( bool bStats )
  {
    m_bEnableStats = bStats;
  }

  QString DataToString();

  Region2D* ObjectFromString(RenderView2D* view, const QString& text);

protected:
  void UpdateWorldCoords();
  int GetRange( double[3][2] );

  vtkSmartPointer<vtkActor2D>   m_actorRect;
  int       m_nX1, m_nX2, m_nY1, m_nY2;
  double    m_dPt[4][3];       // rect in world coordinate

  bool      m_bEnableStats;
};

#endif


