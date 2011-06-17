/**
 * @file  FloatingStatusBar.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/06/17 02:39:27 $
 *    $Revision: 1.7 $
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
#include "FloatingStatusBar.h"
#include "ui_FloatingStatusBar.h"
#include <QDebug>
#include <QDesktopServices>

FloatingStatusBar::FloatingStatusBar(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::FloatingStatusBar)
{
  ui->setupUi(this);
  setWindowFlags(Qt::Tool | Qt::CustomizeWindowHint);
#ifdef Q_WS_MAC
  ui->frame->layout()->setContentsMargins(5, 5, 5, 8);
#elif defined(Q_WS_X11)
  ui->frame->setFrameShape(QFrame::StyledPanel);
#endif
  m_timer = new QTimer( this );
  m_timer->setInterval(250);
  connect( m_timer, SIGNAL(timeout()), this, SLOT(OnProgressTimer()));
}

FloatingStatusBar::~FloatingStatusBar()
{
  delete ui;
}

void FloatingStatusBar::SetProgress( int n )
{
  ui->progressBar->setValue( n );
}

void FloatingStatusBar::ShowProgress()
{
  ui->progressBar->show();
  ui->progressBar->setValue(0);
  this->show();
//  m_timer->start();
  Reposition();
}

void FloatingStatusBar::Reposition()
{
  QWidget* p = parentWidget();
  QSize s = p->size() - this->size();
#ifdef Q_WS_MAC
  this->move(p->geometry().topLeft() + QPoint(0, s.height()) + QPoint(1,-1));
#elif defined(Q_CYGWIN_WIN)
  this->move(p->geometry().topLeft() + QPoint(0, s.height()) + QPoint(0, -6));
#else
  this->move(p->geometry().topLeft() + QPoint(0, s.height()));
#endif
}

void FloatingStatusBar::HideProgress()
{
  m_timer->stop();
  this->hide();
  ui->progressBar->hide();
}

void FloatingStatusBar::OnProgressTimer()
{
  if ( ui->progressBar->value() == 100 )
  {
    ui->progressBar->setValue( 50 );
  }
  else
  {
    ui->progressBar->setValue( ui->progressBar->value()+2 );
  }
}

