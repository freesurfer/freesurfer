/**
 * @brief Coordinate annotation for 2D views.
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

#include "Annotation2D.h"
#include "vtkRenderer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkActor2D.h"
#include "vtkProperty2D.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "vtkPropCollection.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "MyVTKUtils.h"
#include <QDebug>

#define NUMBER_OF_COORD_ANNOTATIONS 6

Annotation2D::Annotation2D( QObject* parent ) : QObject( parent )
{
  m_actorsAll = vtkSmartPointer<vtkPropCollection>::New();
  for ( int i = 0; i < NUMBER_OF_COORD_ANNOTATIONS; i++ )
  {
    m_actorCoordinates[i] = vtkSmartPointer<vtkTextActor>::New();
    m_actorCoordinates[i]->GetPositionCoordinate()->
        SetCoordinateSystemToNormalizedViewport();
    m_actorCoordinates[i]->GetTextProperty()->ShadowOff();
    m_actorCoordinates[i]->GetTextProperty()->SetFontSize(11);
    m_actorCoordinates[i]->GetTextProperty()->ItalicOff();
    //     m_actorCoordinates[i]->GetTextProperty()->SetFontFamilyToTimes();
    //    m_actorCoordinates[i]->SetTextScaleModeToViewport();
    m_actorsAll->AddItem( m_actorCoordinates[i] );
  }
  m_actorCoordinates[0]->SetPosition( 0.01, 0.5 );
  m_actorCoordinates[0]->GetTextProperty()->SetJustificationToLeft();
  m_actorCoordinates[0]->GetTextProperty()->
      SetVerticalJustificationToCentered();
  m_actorCoordinates[1]->SetPosition( 0.5, 0.99 );
  m_actorCoordinates[1]->GetTextProperty()->SetJustificationToCentered();
  m_actorCoordinates[1]->GetTextProperty()->SetVerticalJustificationToTop();
  m_actorCoordinates[2]->SetPosition( 0.99, 0.5 );
  m_actorCoordinates[2]->GetTextProperty()->SetJustificationToRight();
  m_actorCoordinates[2]->GetTextProperty()->
      SetVerticalJustificationToCentered();
  m_actorCoordinates[3]->SetPosition( 0.5, 0.01 );
  m_actorCoordinates[3]->GetTextProperty()->SetJustificationToCentered();
  m_actorCoordinates[3]->GetTextProperty()->SetVerticalJustificationToBottom();

  // indicate slice location
  m_actorCoordinates[4]->SetPosition( 0.99, 0.01 );
  m_actorCoordinates[4]->GetTextProperty()->SetJustificationToRight();
  m_actorCoordinates[4]->GetTextProperty()->SetVerticalJustificationToBottom();

  // indicate slice number
  m_actorCoordinates[5]->SetPosition( 0.99, 0.99 );
  m_actorCoordinates[5]->GetTextProperty()->SetJustificationToRight();
  m_actorCoordinates[5]->GetTextProperty()->SetVerticalJustificationToTop();

  // scale actors
  m_actorScaleLine = vtkSmartPointer<vtkActor2D>::New();
  m_actorScaleLine->SetMapper( vtkSmartPointer<vtkPolyDataMapper2D>::New() );
  int line_w = 1;
#if VTK_MAJOR_VERSION > 7
  line_w = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  m_actorScaleLine->GetProperty()->SetLineWidth( line_w );
  m_actorScaleLine->GetPositionCoordinate()->
      SetCoordinateSystemToNormalizedViewport();
  m_actorScaleLine->SetPosition( 0.05, 0.05 );
  m_actorsAll->AddItem( m_actorScaleLine );

  m_actorScaleTitle = vtkSmartPointer<vtkTextActor>::New();
  m_actorScaleTitle->GetPositionCoordinate()->
      SetCoordinateSystemToNormalizedViewport();
  m_actorScaleTitle->GetTextProperty()->ShadowOff();
  m_actorScaleTitle->GetTextProperty()->SetFontSize(11);
  m_actorScaleTitle->GetTextProperty()->ItalicOff();
  m_actorScaleTitle->GetTextProperty()->SetJustificationToCentered();
  m_actorScaleTitle->GetTextProperty()->SetVerticalJustificationToTop();
  m_actorsAll->AddItem( m_actorScaleTitle );
}

Annotation2D::~Annotation2D()
{}

void Annotation2D::SetAutoScaleText(bool b)
{
  for ( int i = 0; i < NUMBER_OF_COORD_ANNOTATIONS; i++ )
  {
    if (b)
      m_actorCoordinates[i]->SetTextScaleModeToViewport();
    else
      m_actorCoordinates[i]->SetTextScaleModeToNone();
  }
}

void Annotation2D::SetTextSize(int nsize)
{
  for ( int i = 0; i < NUMBER_OF_COORD_ANNOTATIONS; i++ )
  {
    m_actorCoordinates[i]->GetTextProperty()->SetFontSize(nsize);
  }
}

void Annotation2D::Update( vtkRenderer* renderer, int nPlane )
{
  double slicePos[3] = { 0, 0, 0 };
  LayerCollection* lc =
      MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  lc->GetSlicePosition( slicePos );
  bool bHasLayer = ( lc->GetNumberOfLayers() > 0 );

  if ( !bHasLayer )
  {
    return;
  }

  LayerMRI* mri = ( LayerMRI* )lc->GetActiveLayer();
  int nSliceNumber[3];
  double ras[3];
  if ( mri )
  {
    mri->RemapPositionToRealRAS( slicePos, ras );
    mri->RASToOriginalIndex( ras, nSliceNumber );
  }

  double centPos[3];
  double* worigin = lc->GetWorldOrigin();
  double* wsize = lc->GetWorldSize();
  double* wvoxel = lc->GetWorldVoxelSize();
  for ( int i = 0; i < 3; i++ )
  {
    centPos[i] = worigin[i] + wsize[i]/2;
    if ( !mri )
    {
      nSliceNumber[i] = (int)( ( slicePos[i] - worigin[i] ) / wvoxel[i] );
    }
  }

  double pos[4] = { 0, 1, 1, 0 }, tpos = 0;

  QString strg;
  switch ( nPlane )
  {
  case 0:
    renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
    renderer->ViewToWorld( pos[0], pos[1], tpos );
    pos[0] = pos[1];
    pos[1] = tpos;
    renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
    renderer->ViewToWorld( pos[2], pos[3], tpos );
    pos[2] = pos[3];
    pos[3] = tpos;

    tpos = slicePos[0];
    mri->RemapPositionToRealRAS( tpos, pos[0], pos[1],
        slicePos[0], pos[0], pos[1] );
    if ( pos[0] >= 0 )
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("A %.2f", pos[0]).toLatin1().constData() );
    else
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("P %.2f", -pos[0]).toLatin1().constData() );

    if ( pos[1] >= 0 )
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("S %.2f", pos[1]).toLatin1().constData() );
    else
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("I  %.2f", -pos[1]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( tpos, pos[2], pos[3],
        slicePos[0], pos[2], pos[3] );
    if ( pos[2] >= 0 )
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("A %.2f", pos[2]).toLatin1().constData() );
    else
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("P %.2f", -pos[2]).toLatin1().constData() );

    if ( pos[3] >= 0 )
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("S %.2f", pos[3]).toLatin1().constData() );
    else
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("I  %.2f", -pos[3]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( tpos, centPos[1], centPos[2],
        slicePos[0], centPos[1], centPos[2] );
    if ( slicePos[0] >= 0 )
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("R %.2f", slicePos[0]).toLatin1().constData() );
    else
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("L %.2f", -slicePos[0]).toLatin1().constData() );

    break;
  case 1:
    renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
    renderer->ViewToWorld( pos[0], pos[1], tpos );
    pos[1] = tpos;
    renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
    renderer->ViewToWorld( pos[2], pos[3], tpos );
    pos[3] = tpos;

    tpos = slicePos[1];
    mri->RemapPositionToRealRAS( pos[0], tpos, pos[1], pos[0],
        slicePos[1], pos[1] );
    if ( pos[0] >= 0 )
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("R %.2f", pos[0]).toLatin1().constData() );
    else
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("L %.2f", -pos[0]).toLatin1().constData() );

    if ( pos[1] >= 0 )
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("S %.2f", pos[1]).toLatin1().constData() );
    else
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("I  %.2f", -pos[1]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( pos[2], tpos, pos[3], pos[2],
        slicePos[1], pos[3] );

    if ( pos[2] >= 0 )
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("R %.2f", pos[2]).toLatin1().constData() );
    else
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("L %.2f", -pos[2]).toLatin1().constData() );

    if ( pos[3] >= 0 )
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("S %.2f", pos[3]).toLatin1().constData() );
    else
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("I  %.2f", -pos[3]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( centPos[0], tpos, centPos[2], centPos[0],
        slicePos[1], centPos[2] );
    if ( slicePos[1] >= 0 )
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("A %.2f", slicePos[1]).toLatin1().constData() );
    else
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("P %.2f", -slicePos[1]).toLatin1().constData() );

    break;
  case 2:
    renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
    renderer->ViewToWorld( pos[0], pos[1], tpos );
    tpos = 0;
    renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
    renderer->ViewToWorld( pos[2], pos[3], tpos );

    tpos = slicePos[2];
    mri->RemapPositionToRealRAS( pos[0], pos[1], tpos,
        pos[0], pos[1], slicePos[2] );

    if ( pos[0] >= 0 )
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("R %.2f", pos[0]).toLatin1().constData() );
    else
      m_actorCoordinates[0]->SetInput
          ( strg.sprintf("L %.2f", -pos[0]).toLatin1().constData() );

    if ( pos[1] >= 0 )
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("A %.2f", pos[1]).toLatin1().constData() );
    else
      m_actorCoordinates[1]->SetInput
          ( strg.sprintf("P %.2f", -pos[1]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( pos[2], pos[3], tpos,
        pos[2], pos[3], slicePos[2] );

    if ( pos[2] >= 0 )
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("R %.2f", pos[2]).toLatin1().constData() );
    else
      m_actorCoordinates[2]->SetInput
          ( strg.sprintf("L %.2f", -pos[2]).toLatin1().constData() );

    if ( pos[3] >= 0 )
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("A %.2f", pos[3]).toLatin1().constData() );
    else
      m_actorCoordinates[3]->SetInput
          ( strg.sprintf("P %.2f", -pos[3]).toLatin1().constData() );

    mri->RemapPositionToRealRAS( centPos[0], centPos[1], tpos,
        centPos[0], centPos[1], slicePos[2] );
    if ( slicePos[2] >= 0 )
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("S %.2f", slicePos[2]).toLatin1().constData() );
    else
      m_actorCoordinates[4]->SetInput
          ( strg.sprintf("I  %.2f", -slicePos[2]).toLatin1().constData() );

    break;
  }

  // update slice number
  int nOrigPlane = nPlane;
  if ( mri )
  {
    QString ostr = mri->GetOrientationString();
    char ch[3][3] = {"RL", "AP", "IS"};
    for (int i = 0; i < 3; i++)
    {
      if (ostr[i] == ch[nPlane][0] || ostr[i] == ch[nPlane][1])
      {
        nOrigPlane = i;
        break;
      }
    }
  }
  m_actorCoordinates[5]->SetInput
      ( mri ?  QString::number( nSliceNumber[nOrigPlane] ).toLatin1().constData() : "" );

  // update scale line
  double* xy_pos = m_actorScaleLine->GetPosition();
  double w_pos[3], w_pos2[3];
  int nNumOfTicks = 5;
  MyVTKUtils::NormalizedViewportToWorld( renderer,
                                         xy_pos[0],
      xy_pos[1],
      w_pos[0],
      w_pos[1],
      w_pos[2] );
  MyVTKUtils::NormalizedViewportToWorld( renderer,
                                         xy_pos[0] + 0.5,
      xy_pos[1],
      w_pos2[0],
      w_pos2[1],
      w_pos2[2] );
  w_pos[ nPlane ] = w_pos2[ nPlane ] = 0;
  double d = 0.5 /
      sqrt( vtkMath::Distance2BetweenPoints( w_pos, w_pos2 ) ) * 10;
  QString title = "1 cm";
  if ( d >= 0.5 - xy_pos[0] )
  {
    d /= 2;
    title = "5 mm";
    if ( d >= 0.5 - xy_pos[0] )
    {
      d *= 0.4;
      title = "2 mm";
      nNumOfTicks = 2;
      if ( d >= 0.5 - xy_pos[0] )
      {
        d /= 2;
        title = "1 mm";
        nNumOfTicks = 5;
        if ( d >= 0.5 - xy_pos[0] )
        {
          d /= 2;
          title = "500 um";
          nNumOfTicks = 5;
          if ( d >= 0.5 - xy_pos[0] )
          {
            d *= 0.4;
            title = "200 um";
            nNumOfTicks = 2;
            if ( d >= 0.5 - xy_pos[0] )
            {
              d /= 2;
              title = "100 um";
              nNumOfTicks = 5;
              if ( d >= 0.5 - xy_pos[0] )
              {
                d /= 2;
                title = "50 um";
                nNumOfTicks = 5;
                if ( d >= 0.5 - xy_pos[0] )
                {
                  d *= 0.4;
                  title = "20 um";
                  nNumOfTicks = 2;
                  if ( d >= 0.5 - xy_pos[0] )
                  {
                    d *= 0.5;
                    title = "10 um";
                    nNumOfTicks = 5;
                    if ( d >= 0.5 - xy_pos[0] )
                    {
                      d *= 0.5;
                      title = "5 um";
                      nNumOfTicks = 5;
                      if ( d >= 0.5 - xy_pos[0] )
                      {
                        d *= 0.4;
                        title = "2 um";
                        nNumOfTicks = 2;
                        if ( d >= 0.5 - xy_pos[0] )
                        {
                          d *= 0.5;
                          title = "1 um";
                          nNumOfTicks = 5;
                          if ( d >= 0.5 - xy_pos[0] )
                          {
                            d *= 0.5;
                            title = "500 nm";
                            nNumOfTicks = 5;
                            if ( d >= 0.5 - xy_pos[0] )
                            {
                              d *= 0.4;
                              title = "200 nm";
                              nNumOfTicks = 2;
                              if ( d >= 0.5 - xy_pos[0] )
                              {
                                d *= 0.5;
                                title = "100 nm";
                                nNumOfTicks = 5;
                                if ( d >= 0.5 - xy_pos[0] )
                                {
                                  d *= 0.5;
                                  title = "50 nm";
                                  nNumOfTicks = 5;
                                  if ( d >= 0.5 - xy_pos[0] )
                                  {
                                    d *= 0.4;
                                    title = "20 nm";
                                    nNumOfTicks = 2;
                                    if ( d >= 0.5 - xy_pos[0] )
                                    {
                                      d *= 0.5;
                                      title = "10 nm";
                                      nNumOfTicks = 1;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  UpdateScaleActors( d, nNumOfTicks,  title.toLatin1().constData() );
}

void Annotation2D::UpdateScaleActors( double length,
                                      int nNumOfTicks,
                                      const char* title )
{
  // scale line
  double* pos = m_actorScaleLine->GetPosition();
  double tick_len = 0.007;
  vtkPoints* Pts = vtkPoints::New();
  Pts->InsertNextPoint( 0, 0, 0 );
  Pts->InsertNextPoint( length, 0, 0 );
  vtkCellArray* Lines = vtkCellArray::New();
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 0 );
  Lines->InsertCellPoint( 1 );
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 1 );
  Lines->InsertCellPoint( 0 );

  int n = 2;
  double step = length / nNumOfTicks;
  for ( int i = 0; i <= nNumOfTicks; i++ )
  {
    Pts->InsertNextPoint( step*i, 0, 0 );
    Pts->InsertNextPoint( step*i, -tick_len, 0 );
    Lines->InsertNextCell( 2 );
    Lines->InsertCellPoint( n++ );
    Lines->InsertCellPoint( n++ );
  }

  vtkPolyData* poly = vtkPolyData::New();
  poly->SetPoints(Pts);
  poly->SetLines(Lines);
  Pts->Delete();
  Lines->Delete();

  vtkCoordinate* normCoords = vtkCoordinate::New();
  normCoords->SetCoordinateSystemToNormalizedViewport();

  vtkPolyDataMapper2D* pMapper = vtkPolyDataMapper2D::New();
#if VTK_MAJOR_VERSION > 5
  pMapper->SetInputData( poly );
#else
  pMapper->SetInput(poly);
#endif
  pMapper->SetTransformCoordinate(normCoords);
  poly->Delete();
  normCoords->Delete();

  m_actorScaleLine->SetMapper(pMapper);
  pMapper->Delete();

  // title
  m_actorScaleTitle->SetPosition( pos[0] + length / 2, pos[1] - 0.01 );
  m_actorScaleTitle->SetInput( title );
}

void Annotation2D::AppendAnnotations( vtkRenderer* renderer, bool bScaleBar )
{
  if (!bScaleBar)
  {
    for ( int i = 0; i < NUMBER_OF_COORD_ANNOTATIONS; i++ )
    {
      renderer->AddViewProp( m_actorCoordinates[i] );
    }
  }
  else
  {
    renderer->AddViewProp( m_actorScaleLine );
    renderer->AddViewProp( m_actorScaleTitle );
  }
}

void Annotation2D::ShowScaleLine( bool bShow )
{
  m_actorScaleLine->SetVisibility( bShow?1:0 );
  m_actorScaleTitle->SetVisibility( bShow?1:0 );
}

bool Annotation2D::GetShowScaleLine()
{
  return m_actorScaleLine->GetVisibility() > 0;
}

void Annotation2D::Show( bool bShow )
{
  vtkProp* prop = NULL;
  m_actorsAll->InitTraversal();
  while ( ( prop = m_actorsAll->GetNextProp() ) )
  {
    prop->SetVisibility( bShow?1:0 );
  }
}

bool Annotation2D::IsVisible()
{
  m_actorsAll->InitTraversal();
  vtkProp* prop = m_actorsAll->GetNextProp();
  if ( prop )
  {
    return (prop->GetVisibility() > 0);
  }
  else
  {
    return false;
  }
}

QColor Annotation2D::GetColor()
{
  double c[3];
  m_actorScaleLine->GetProperty()->GetColor(c);
  return QColor::fromRgbF(c[0], c[1], c[2]);
}

void Annotation2D::SetColor(const QColor &c)
{
  vtkProp* prop = NULL;
  m_actorsAll->InitTraversal();
  while ( ( prop = m_actorsAll->GetNextProp() ) )
  {
    vtkActor2D* actor = vtkActor2D::SafeDownCast(prop);
    if (actor)
      actor->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  }
  emit Updated();
}
