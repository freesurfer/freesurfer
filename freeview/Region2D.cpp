/**
 * @file  Region2D.cpp
 * @brief Region2D.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2016/03/17 16:25:26 $
 *    $Revision: 1.12 $
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
  m_actorText->SetTextScaleModeToNone();
  m_actorText->GetTextProperty()->SetColor( 1, 1, 1 );
  m_actorText->GetTextProperty()->ShadowOn();
  m_actorText->GetTextProperty()->SetFontSize( 20 );
  m_actorText->GetTextProperty()->SetFontFamilyToTimes();
  m_actorText->GetTextProperty()->SetJustificationToCentered();
}

Region2D::~Region2D()
{}

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

