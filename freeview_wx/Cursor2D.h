/**
 * @file  Cursor2D.h
 * @brief Cursor for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
 *    $Revision: 1.1 $
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

#ifndef Cursor2D_h
#define Cursor2D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"
#include <wx/colour.h>
#include <vector>

class vtkActor2D;
class vtkRenderer;
class vtkActor;
class vtkCursor2D;
class RenderView2D;

class Cursor2D
{
public:
  Cursor2D( RenderView2D* view );
  virtual ~Cursor2D();

  enum CursorStyle { CS_Short = 0, CS_Long };
  
  void SetPosition( double* pos, bool bConnectPrevious = false );
  void SetPosition2( double* pos);

  double* GetPosition();
  void GetPosition( double* pos );

  void SetInterpolationPoints( std::vector<double> pts );

  std::vector<double> GetInterpolationPoints()
  {
    return m_dInterpolationPoints;
  }

  void ClearInterpolationPoints()
  {
    m_dInterpolationPoints.clear();
  }

  void GetColor( double* rgb );
  void SetColor( double r, double g, double b );

  wxColour GetColor();
  void SetColor( const wxColour& color );

  int GetRadius();
  void SetRadius( int nPixels );

  void Update( bool bConnectPrevious = false );

  void AppendActor( vtkRenderer* renderer );

  void Show( bool bShow = true );

  bool IsShown();
  
  int GetStyle()
  {
    return m_nStyle;
  }
  
  void SetStyle( int nStyle );

private:
  vtkSmartPointer<vtkActor2D> m_actorCursor;

  RenderView2D* m_view;

  double  m_dPosition[3];
  double  m_dPosition2[3];
  std::vector<double> m_dInterpolationPoints;

  int   m_nRadius;
  int   m_nStyle;
};

#endif


