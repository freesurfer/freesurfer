/**
 * @brief Region2D.
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

#include "Region2D.h"
#include "RenderView2D.h"
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <QClipboard>
#include <QApplication>

Region2D::Region2D( RenderView2D* view ) :
  QObject( view ),
  m_view( view )
{
  m_actorText = vtkSmartPointer<vtkTextActor>::New();
  m_actorText->GetTextProperty()->SetColor( 1, 1, 1 );
  m_actorText->GetTextProperty()->ShadowOn();
  m_actorText->GetTextProperty()->SetFontSize( 20 );
  m_actorText->GetTextProperty()->SetFontFamilyToTimes();
  m_actorText->GetTextProperty()->SetJustificationToCentered();
  SetTextSize(view->GetTextSize());
  SetAutoScaleText(view->GetAutoScaleText());
}

Region2D::~Region2D()
{}

vtkTextActor* Region2D::GetTextActor()
{
  return m_actorText;
}

void Region2D::SetTextSize(int nsize)
{
  m_actorText->GetTextProperty()->SetFontSize( nsize*1.5 );
}

void Region2D::SetAutoScaleText(bool b)
{
  if (b)
    m_actorText->SetTextScaleModeToViewport();
  else
    m_actorText->SetTextScaleModeToNone();
}

void Region2D::UpdateStats()
{
  //  QString text = m_actorText->GetInput();
  //  QClipboard *clipboard = QApplication::clipboard();
  //  clipboard->setText(text);

  emit StatsUpdated();
}

void Region2D::UpdateSlicePosition( int nPlane, double pos )
{
  Q_UNUSED(nPlane);
  Q_UNUSED(pos);
  Update();
}

void Region2D::Show( bool bShow )
{
  Q_UNUSED(bShow);
}

