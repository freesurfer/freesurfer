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
#include "DialogSavePointSet.h"
#include "ui_DialogSavePointSet.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
#include "LayerPropertyPointSet.h"

DialogSavePointSet::DialogSavePointSet(QWidget *parent) :
  QDialog(parent),
  m_bRemind(false),
  ui(new Ui::DialogSavePointSet)
{
  ui->setupUi(this);
}

DialogSavePointSet::~DialogSavePointSet()
{
  delete ui;
}

void DialogSavePointSet::SetFileName(const QString &fn_in, int type)
{
  if (fn_in.isEmpty())
    return;

  QString fn = fn_in;
  QString suffix = QFileInfo(fn).suffix();
  if (suffix == "label" || suffix == "dat" || suffix == "json")
    fn = fn.left(fn.length()-suffix.length()-1);
  if (type == LayerPropertyPointSet::WayPoint)
    fn += ".label";
  else if (type == LayerPropertyPointSet::ControlPoint)
    fn += ".dat";
  else
    fn += ".json";
  ui->lineEditFileName->setText(fn);
  ui->lineEditFileName->setCursorPosition(fn.length());
}

QString DialogSavePointSet::GetFileName()
{
  return ui->lineEditFileName->text().trimmed();
}

void DialogSavePointSet::OnOK()
{
  QString fn = GetFileName();
  if (fn.isEmpty())
  {
    QMessageBox::information(this, "Error", "Please enter a file name to save.");
    return;
  }
  else
  {
    if (QFile(fn).exists())
    {
      if (QMessageBox::question(this, "Confirm",
                                tr("%1 already exists. Do you want to replace it?").arg(QFileInfo(fn).fileName()),
                                QMessageBox::No, QMessageBox::Yes) == QMessageBox::No)
        return;
    }
    if (m_bRemind && GetType() != LayerPropertyPointSet::Enhanced)
    {
      if (QMessageBox::question(this, "Warning", "Point set data contains enhanced information such as comments. "
                                "Saving it in legacy formats will lose all the enhanced info. Do you still want to continue?",
                                QMessageBox::No, QMessageBox::Yes)
          == QMessageBox::No)
        return;
    }
    accept();
  }
}

void DialogSavePointSet::OnOpen()
{
  QString old_fn = GetFileName();
  QString fn = QFileDialog::getSaveFileName(this, "Select File To Save",
                                            (old_fn.isEmpty()?m_strLastDir:old_fn),
                                            "All Files (*)", NULL, QFileDialog::DontConfirmOverwrite);
  if (!fn.isEmpty())
  {
    SetFileName(fn, GetType());
  }
}

void DialogSavePointSet::SetType(int nType)
{
  ui->radioButtonControlPoint->setChecked(nType == LayerPropertyPointSet::ControlPoint);
  ui->radioButtonWayPoint->setChecked(nType == LayerPropertyPointSet::WayPoint);
  ui->radioButtonEnhanced->setChecked(nType == LayerPropertyPointSet::Enhanced);
  if (nType == LayerPropertyPointSet::Enhanced)
    m_bRemind = true;
}

int DialogSavePointSet::GetType()
{
  if ( ui->radioButtonControlPoint->isChecked())
  {
    return LayerPropertyPointSet::ControlPoint;
  }
  else if (ui->radioButtonWayPoint->isChecked())
  {
    return LayerPropertyPointSet::WayPoint;
  }
  else
    return LayerPropertyPointSet::Enhanced;
}

void DialogSavePointSet::OnRadioButtonFileType()
{
  SetFileName(ui->lineEditFileName->text().trimmed(), GetType());
}
