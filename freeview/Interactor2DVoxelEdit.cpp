/**
 * @brief Interactor for editing voxel in 2D render view.
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

#include "Interactor2DVoxelEdit.h"
#include <QDebug>
#include <QMouseEvent>
#include <QKeyEvent>

Interactor2DVoxelEdit::Interactor2DVoxelEdit( QObject* parent ) :
  Interactor2DVolumeEdit( "MRI", parent )
{}

bool Interactor2DVoxelEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* view )
{
  return Interactor2DVolumeEdit::ProcessMouseDownEvent(event, view);
}

bool Interactor2DVoxelEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* view )
{
  return Interactor2DVolumeEdit::ProcessMouseUpEvent(event, view);
}

bool Interactor2DVoxelEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  return Interactor2DVolumeEdit::ProcessKeyDownEvent(event, renderview);
}

bool Interactor2DVoxelEdit::ProcessKeyUpEvent( QKeyEvent* event, RenderView* renderview )
{
  return Interactor2DVolumeEdit::ProcessKeyUpEvent(event, renderview);
}
