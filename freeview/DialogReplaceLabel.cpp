/**
 * @brief Dialog window to replace label.
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

#include "DialogReplaceLabel.h"
#include "ui_DialogReplaceLabel.h"
#include <QMessageBox>

DialogReplaceLabel::DialogReplaceLabel(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogReplaceLabel)
{
  ui->setupUi(this);
}

DialogReplaceLabel::~DialogReplaceLabel()
{
  delete ui;
}

void DialogReplaceLabel::OnReplace()
{
  bool bOK;
  ui->lineEditOriginalValue->text().toDouble(&bOK);
  if (!bOK)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid original value.");
    return;
  }
  ui->lineEditNewValue->text().toDouble(&bOK);
  if (!bOK)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid new value.");
    return;
  }
  this->accept();
}

double DialogReplaceLabel::GetOriginalValue()
{
  return ui->lineEditOriginalValue->text().toDouble();
}

double DialogReplaceLabel::GetNewValue()
{
  return ui->lineEditNewValue->text().toDouble();
}

bool DialogReplaceLabel::ReplaceSingleSlice()
{
  return ui->radioButtonCurrentSlice->isChecked();
}
