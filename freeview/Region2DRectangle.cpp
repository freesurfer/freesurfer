/**
 * @file  Region2DRectangle.cpp
 * @brief Region2DRectangle.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.15 $
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

#include "Region2DRectangle.h"
#include "RenderView2D.h"
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCoordinate.h>
#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include "LayerMRI.h"
#include <QDebug>

Region2DRectangle::Region2DRectangle( RenderView2D* view ) :
  Region2D( view )
{
  m_bEnableStats = true;

  m_actorRect = vtkSmartPointer<vtkActor2D>::New();
  m_actorRect->GetProperty()->SetColor( 1, 0, 0 );
  m_actorRect->GetProperty()->SetOpacity( 0.6 );
//  m_actorRect->VisibilityOff();

  Highlight( true );
}

Region2DRectangle::~Region2DRectangle()
{}


void Region2DRectangle::SetRect( int left, int top, int w, int h )
{
  m_nX1 = left;
  m_nY1 = top;
  m_nX2 = left+w;
  m_nY2 = top+h;

  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::UpdateWorldCoords()
{
  m_view->MousePositionToRAS( m_nX1, m_nY1, m_dPt[0] );
  m_view->MousePositionToRAS( m_nX1, m_nY2, m_dPt[1] );
  m_view->MousePositionToRAS( m_nX2, m_nY2, m_dPt[2] );
  m_view->MousePositionToRAS( m_nX2, m_nY1, m_dPt[3] );
}

void Region2DRectangle::SetTopLeft( int left, int top )
{
  m_nX1 = left;
  m_nY1 = top;

  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::SetBottomRight( int right, int bottom )
{
  m_nX2 = right;
  m_nY2 = bottom;

  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::Update()
{
  if ( !m_actorRect->GetVisibility() )
  {
    return;
  }

  double pt0[3], pt1[3], pt2[3], pt3[3];

  m_view->WorldToViewport( m_dPt[0][0], m_dPt[0][1], m_dPt[0][2], pt0[0], pt0[1], pt0[2] );
  m_view->WorldToViewport( m_dPt[1][0], m_dPt[1][1], m_dPt[1][2], pt1[0], pt1[1], pt1[2] );
  m_view->WorldToViewport( m_dPt[2][0], m_dPt[2][1], m_dPt[2][2], pt2[0], pt2[1], pt2[2] );
  m_view->WorldToViewport( m_dPt[3][0], m_dPt[3][1], m_dPt[3][2], pt3[0], pt3[1], pt3[2] );

  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> poly = vtkSmartPointer<vtkCellArray>::New();
  pts->InsertNextPoint( pt0 );
  pts->InsertNextPoint( pt1 );
  pts->InsertNextPoint( pt2 );
  pts->InsertNextPoint( pt3 );
  vtkIdType face[4] = { 0, 1, 2, 3 };
  poly->InsertNextCell( 4, face );

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( pts );
  polydata->SetPolys( poly );
  vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
  mapper->SetInput( polydata );

  vtkSmartPointer<vtkCoordinate> coords = vtkSmartPointer<vtkCoordinate>::New();
  coords->SetCoordinateSystemToViewport();
  mapper->SetTransformCoordinate( coords );
  m_actorRect->SetMapper( mapper );

  if ( m_bEnableStats )
  {
    UpdateStats();

    double pt[3];
    for ( int i = 0; i < 3; i++ )
    {
      pt[i] = ( pt0[i] + pt2[i] ) / 2;
    }
    m_actorText->SetInput( m_strShortStats.toAscii().constData() );
    m_actorText->SetPosition( pt );
  }
}

void Region2DRectangle::UpdateStats()
{
  LayerMRI* layer = m_view->GetFirstNonLabelVolume();
  if ( layer )
  {
    double mean, sd;
    int count;
    layer->GetVoxelStatsRectangle( m_dPt[0], m_dPt[2], m_view->GetViewPlane(), &mean, &sd, &count );
    char ch[1000];
    sprintf( ch, "%.2f +/-%.2f", mean, sd );
    m_actorText->SetInput( ch );
    m_strShortStats = ch;

    m_strsLongStats.clear();
    m_strsLongStats.push_back( QString("Number of voxels: %1").arg(count) );
    m_strsLongStats.push_back( QString("Mean: ") + m_strShortStats );
  }

  Region2D::UpdateStats();
}

void Region2DRectangle::AppendProp( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorRect );
  renderer->AddViewProp( m_actorText );
}

void Region2DRectangle::Show( bool bShow )
{
  m_actorRect->SetVisibility( bShow?1:0 );
  Update();
}

void Region2DRectangle::GetWorldPoint( int nIndex, double* pt )
{
  if ( nIndex >= 0 && nIndex <= 3 )
  {
    pt[0] = m_dPt[nIndex][0];
    pt[1] = m_dPt[nIndex][1];
    pt[2] = m_dPt[nIndex][2];
  }
}

bool Region2DRectangle::Contains( int nX, int nY, int* nIndexOut )
{
  // first calculate the threshold distance in world space
  double pt1[3], pt2[3];
  m_view->MousePositionToRAS( 0, 0, pt1 );
  m_view->MousePositionToRAS( 5, 5, pt2 );
  double dTh2 = vtkMath::Distance2BetweenPoints( pt1, pt2 );

  // calculate the hit point in world space
  double pt[3];
  m_view->MousePositionToRAS( nX, nY, pt );
  int nPlane = m_view->GetViewPlane();
  pt[nPlane] = m_dPt[0][nPlane];

  for ( int i = 0; i < 4; i++ )
  {
    if ( vtkMath::Distance2BetweenPoints( pt, m_dPt[i] ) < dTh2 )
    {
      if ( nIndexOut )
      {
        *nIndexOut = i;
      }
      return true;
    }
  }

  if ( nIndexOut )
  {
    *nIndexOut = -1;
  }

  double range[3][2];
  GetRange( range );
  if ( pt[0] >= range[0][0] && pt[0] <= range[0][1] &&
       pt[1] >= range[1][0] && pt[1] <= range[1][1] &&
       pt[2] >= range[2][0] && pt[2] <= range[2][1] )
  {
    return true;
  }
  else
  {
    return false;
  }
}

// return the index of the screen depth dimension
int Region2DRectangle::GetRange( double range[3][2] )
{
  for ( int i = 0; i < 3; i++ )
  {
    range[i][0] = range[i][1] = m_dPt[0][i];
  }

  for ( int i = 0; i < 4; i++ )
  {
    for ( int j = 0; j < 3; j++ )
    {
      if ( range[j][0] > m_dPt[i][j] )
      {
        range[j][0] = m_dPt[i][j];
      }
      else if ( range[j][1] < m_dPt[i][j] )
      {
        range[j][1] = m_dPt[i][j];
      }
    }
  }

  for ( int i = 0; i < 3; i++ )
  {
    if ( fabs( range[i][0] - range[i][1] ) < 0.000001 )
    {
      return i;
    }
  }

  return 0;
}

void Region2DRectangle::Highlight( bool bHighlight )
{
  m_actorRect->GetProperty()->SetColor( bHighlight?0:1, 1, 0);
}

void Region2DRectangle::Offset( int x, int y )
{
  double pt0[3], pt[3];
  m_view->MousePositionToRAS( 0, 0, pt0 );
  m_view->MousePositionToRAS( x, y, pt );
  for ( int j = 0; j < 4; j++ )
  {
    for ( int i = 0; i < 3; i++ )
    {
      m_dPt[j][i] += ( pt[i] - pt0[i] );
    }
  }
  Update();
}

void Region2DRectangle::UpdatePoint( int nIndex, int x, int y )
{
  /*
  if ( nIndex >= 0 && nIndex < 4 )
  {
   m_view->MousePositionToRAS( nX, nY, m_dPt[nIndex] );

   Update();
  }
  */
  if ( nIndex >= 0 && nIndex < 4 )
  {
    int nDiag = (nIndex+2)%4;
    int xp, yp;
    m_view->WorldToScreen( m_dPt[nDiag][0], m_dPt[nDiag][1], m_dPt[nDiag][2], xp, yp );
    int nX1, nY1, nX2, nY2;
    if ( nIndex == 0 )
    {
      nX1 = x;
      nY1 = y;
      nX2 = xp;
      nY2 = yp;
    }
    else if ( nIndex == 1 )
    {
      nX1 = x;
      nY1 = yp;
      nX2 = xp;
      nY2 = y;
    }
    else if ( nIndex == 2 )
    {
      nX1 = xp;
      nY1 = yp;
      nX2 = x;
      nY2 = y;
    }
    else
    {
      nX1 = xp;
      nY1 = y;
      nX2 = x;
      nY2 = yp;
    }
    SetRect( nX1, nY1, nX2-nX1, nY2-nY1 );
  }
}

void Region2DRectangle::UpdateSlicePosition( int nPlane, double pos )
{
  for ( int i = 0; i < 4; i++ )
  {
    m_dPt[i][nPlane] = pos;
  }

  Region2D::UpdateSlicePosition( nPlane, pos );
}

QString Region2DRectangle::DataToString()
{
  QString strg = "FreeView:Region2DRectangle:";
  for (int i = 0; i < 4; i++)
    strg += QString("%1,%2,%3,").arg(m_dPt[i][0]).arg(m_dPt[i][1]).arg(m_dPt[i][2]);
  strg.chop(1);
  return strg;
}

Region2D* Region2DRectangle::ObjectFromString(RenderView2D* view, const QString& text)
{
  QString head = "FreeView:Region2DRectangle:";
  if (text.indexOf(head) != 0)
    return NULL;

  QStringList list = text.mid(head.size()).split(",");
  if (list.size() < 12)
    return NULL;

  double dval[12];
  bool bOK = true;
  dval[0] = list[0].toDouble(&bOK);
  int i = 1;
  while( bOK && i < 12 )
  {
    dval[i] = list[i].toDouble(&bOK);
    i++;
  }
  Region2DRectangle* reg = NULL;
  if (bOK)
  {
    reg = new Region2DRectangle(view);
    for (int i = 0; i < 12; i++)
      reg->m_dPt[i/3][i%3] = dval[i];
    reg->Update();
  }
  return reg;
}
