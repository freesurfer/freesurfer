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
#include "DialogGradientFilter.h"
#include "ui_DialogGradientFilter.h"

DialogGradientFilter::DialogGradientFilter(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogGradientFilter)
{
  ui->setupUi(this);
}

DialogGradientFilter::~DialogGradientFilter()
{
  delete ui;
}

void DialogGradientFilter::SetSmoothing(bool smooth)
{
  ui->checkBoxSmooth->setChecked(smooth);
}

bool DialogGradientFilter::GetSmoothing()
{
  return ui->checkBoxSmooth->isChecked();
}

void DialogGradientFilter::SetSD(double val)
{
  ui->doubleSpinBoxSD->setValue(val);
}

double DialogGradientFilter::GetSD()
{
  return ui->doubleSpinBoxSD->value();
}
