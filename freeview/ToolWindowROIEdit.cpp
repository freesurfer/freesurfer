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
 */
#include "ToolWindowROIEdit.h"
#include "ui_ToolWindowROIEdit.h"
#include "Interactor2DROIEdit.h"
#include "RenderView2D.h"
#include "MainWindow.h"
#include "BrushProperty.h"
#include <QSettings>
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif

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
  MainWindow* wnd = MainWindow::GetMainWindow();
  connect( ui->spinBoxBrushSize, SIGNAL(valueChanged(int)), wnd->GetBrushProperty(), SLOT(SetBrushSize(int)));

  UpdateWidgets();

#ifdef Q_OS_MAC
  if (MacHelper::IsDarkMode())
  {
      ui->actionFreeHand->setIcon(MacHelper::InvertIcon(ui->actionFreeHand->icon(), QSize(), true));
      ui->actionPolyLine->setIcon(MacHelper::InvertIcon(ui->actionPolyLine->icon(), QSize(), true));
  }
#endif
}

ToolWindowROIEdit::~ToolWindowROIEdit()
{
  QSettings settings;
  settings.setValue("ToolWindowROIEdit/Position", pos()-this->parentWidget()->pos());

  delete ui;
}

void ToolWindowROIEdit::showEvent(QShowEvent* event)
{
  Q_UNUSED(event);
  static bool bFirstTime = true;
  if ( bFirstTime )
  {
    this->move( parentWidget()->pos() + QPoint(20,100) );
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

  MainWindow* wnd = MainWindow::GetMainWindow();
  RenderView2D* view = (RenderView2D*)wnd->GetRenderView( 0 );
  ui->actionFill->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Fill );
  ui->actionLiveWire->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Livewire );
  ui->actionFreeHand->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Freehand );
  ui->actionPolyLine->setChecked( view->GetAction() == Interactor2DROIEdit::EM_Polyline );
  ui->spinBoxBrushSize->setValue( wnd->GetBrushProperty()->GetBrushSize() );

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
