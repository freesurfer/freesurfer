/**
 * @file  ToolWindowROIEdit.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/04/11 19:46:21 $
 *    $Revision: 1.4.2.3 $
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
#include "ToolWindowROIEdit.h"
#include "ui_ToolWindowROIEdit.h"
#include "Interactor2DROIEdit.h"
#include "RenderView2D.h"
#include "MainWindow.h"
#include <QSettings>

ToolWindowROIEdit::ToolWindowROIEdit(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::ToolWindowROIEdit)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction( ui->actionFill );
  ag->addAction( ui->actionFreeHand );
  ag->addAction( ui->actionLiveWire );
  ag->addAction( ui->actionPolyLine );
  ag->setExclusive( true );
  ui->actionFill->setData( Interactor2DROIEdit::EM_Fill );
  ui->actionFreeHand->setData( Interactor2DROIEdit::EM_Freehand );
  ui->actionLiveWire->setData( Interactor2DROIEdit::EM_Livewire );
  ui->actionPolyLine->setData( Interactor2DROIEdit::EM_Polyline );
  connect( ag, SIGNAL(triggered(QAction*)), this, SLOT(OnEditMode(QAction*)) );

  UpdateWidgets();
}

ToolWindowROIEdit::~ToolWindowROIEdit()
{
  QSettings settings;
  settings.setValue("ToolWindowROIEdit/Position", pos()-this->parentWidget()->pos());

  delete ui;
}

void ToolWindowROIEdit::showEvent(QShowEvent* event)
{
  static bool bFirstTime = true;
  if ( bFirstTime )
  {
    QSettings settings;
    QVariant v = settings.value( "ToolWindowROIEdit/Position", QPoint( 200, 20 ) );
    this->move( parentWidget()->pos() + v.toPoint() );
    bFirstTime = false;
  }
}

void ToolWindowROIEdit::UpdateWidgets( )
{
  QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( true );
  }

  RenderView2D* view = (RenderView2D*)MainWindow::GetMainWindow()->GetRenderView( 0 );
  ui->actionFill->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Fill );
  ui->actionLiveWire->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Livewire );
  ui->actionFreeHand->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Freehand );
  ui->actionPolyLine->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Polyline );

  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( false );
  }
}

void ToolWindowROIEdit::OnEditMode(QAction *act)
{
  MainWindow::GetMainWindow()->SetAction( act->data().toInt() );
  UpdateWidgets();
}
