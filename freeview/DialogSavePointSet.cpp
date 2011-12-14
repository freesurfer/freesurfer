/**
 * @file  DialogSavePointSet.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/14 17:13:44 $
 *    $Revision: 1.6 $
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
#include "DialogSavePointSet.h"
#include "ui_DialogSavePointSet.h"
#include <QMessageBox>
#include <QFileDialog>
#include "LayerPropertyPointSet.h"

DialogSavePointSet::DialogSavePointSet(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSavePointSet)
{
  ui->setupUi(this);
}

DialogSavePointSet::~DialogSavePointSet()
{
  delete ui;
}

void DialogSavePointSet::SetFileName(const QString &fn)
{
  ui->lineEditFileName->setText(fn);
  ui->lineEditFileName->setCursorPosition(fn.length());
}

QString DialogSavePointSet::GetFileName()
{
  return ui->lineEditFileName->text().trimmed();
}

void DialogSavePointSet::OnOK()
{
  if (GetFileName().isEmpty())
  {
    QMessageBox::information(this, "Error", "Please enter a file name to save.");
    return;
  }
  else
  {
    accept();
  }
}

void DialogSavePointSet::OnOpen()
{
  QString old_fn = GetFileName();
  QString fn = QFileDialog::getSaveFileName(this, "Select File To Save",
               (old_fn.isEmpty()?m_strLastDir:old_fn),                                              "All Files (*)");
  if (!fn.isEmpty())
  {
    SetFileName(fn);
  }
}

void DialogSavePointSet::SetType(int nType)
{
  ui->radioButtonControlPoint->setChecked(nType == LayerPropertyPointSet::ControlPoint);
  ui->radioButtonWayPoint->setChecked(nType == LayerPropertyPointSet::WayPoint);
}

int DialogSavePointSet::GetType()
{
  if ( ui->radioButtonControlPoint->isChecked())
  {
    return LayerPropertyPointSet::ControlPoint;
  }
  else
  {
    return LayerPropertyPointSet::WayPoint;
  }
}
