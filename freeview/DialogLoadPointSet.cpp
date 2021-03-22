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
#include "DialogLoadPointSet.h"
#include "ui_DialogLoadPointSet.h"
#include "LayerPropertyPointSet.h"
#include "MyUtils.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>

DialogLoadPointSet::DialogLoadPointSet(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadPointSet)
{
  ui->setupUi(this);
}

DialogLoadPointSet::~DialogLoadPointSet()
{
  delete ui;
}

void DialogLoadPointSet::OnOK()
{
  if ( GetFileNames().isEmpty() )
  {
    QMessageBox::warning(
          this, "Error",
          "Point set file name can not be empty." );
    return;
  }
  accept();
}

void DialogLoadPointSet::OnButtonOpen()
{
  QStringList fns = QFileDialog::getOpenFileNames(
        this,
        "Select point set files",
        m_strLastDir,
        "All files (*)" );
  if ( !fns.isEmpty())
  {
    m_strLastDir = QFileInfo(fns[0]).canonicalPath();
    for (int i = 0; i < fns.size(); i++)
    {
      fns[i] = MyUtils::Win32PathProof(fns[i]);
    }
    ui->lineEditFileName->setText( fns.join(";") );
    ui->lineEditFileName->setCursorPosition(ui->lineEditFileName->text().size());
  }
}

int DialogLoadPointSet::GetPointSetType()
{
  if (ui->radioButtonControlPoint->isChecked() )
  {
    return LayerPropertyPointSet::ControlPoint;
  }
  else if ( ui->radioButtonWayPoint->isChecked() )
  {
    return LayerPropertyPointSet::WayPoint;
  }
  else
  {
    return -1;
  }
}

QStringList DialogLoadPointSet::GetFileNames()
{
  QStringList fns = ui->lineEditFileName->text().trimmed().split(";", QString::SkipEmptyParts);
  for (int i = 0; i < fns.size(); i++)
  {
    fns[i] = MyUtils::CygwinPathProof(fns[i]);
  }
  return fns;
}
