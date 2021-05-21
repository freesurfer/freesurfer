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
#include "DialogVolumeFilter.h"
#include "ui_DialogVolumeFilter.h"
#include "VolumeFilter.h"
#include <QMessageBox>

DialogVolumeFilter::DialogVolumeFilter(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogVolumeFilter)
{
  ui->setupUi(this);
}

DialogVolumeFilter::~DialogVolumeFilter()
{
  delete ui;
}

void DialogVolumeFilter::SetFilter( VolumeFilter* filter )
{
  m_filter = filter;
  if ( filter )
  {
    ui->spinBoxKernelSize->setValue( filter->GetKernelSize() );
  }
}

int DialogVolumeFilter::GetKernelSize()
{
  return ui->spinBoxKernelSize->value();
}

void DialogVolumeFilter::SetSigma( double dvalue )
{
  ui->lineEditSigma->setText( QString::number(dvalue ) );
}

double DialogVolumeFilter::GetSigma()
{
  double dvalue = 0;
  dvalue = ui->lineEditSigma->text().toDouble();
  return dvalue;
}

void DialogVolumeFilter::ShowSigma( bool bShow )
{
  ui->labelSigma->setVisible(bShow);
  ui->lineEditSigma->setVisible(bShow);
}

void DialogVolumeFilter::OnOK()
{
  if ( GetKernelSize() <= 0 )
  {
    QMessageBox::warning( this, "Error", "Kernel size must be greater than 0.");
    return;
  }
  else if ( ui->labelSigma->isVisible() && GetSigma() <= 0 )
  {
    QMessageBox::warning(this, "Error", "Sigma must be greater than 0.");
    return;
  }
  accept();
}
