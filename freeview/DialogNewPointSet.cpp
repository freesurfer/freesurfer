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
#include "DialogNewPointSet.h"
#include "ui_DialogNewPointSet.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "LayerPropertyPointSet.h"
#include <QMessageBox>

DialogNewPointSet::DialogNewPointSet(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogNewPointSet)
{
  ui->setupUi(this);
  LayerCollection* col_mri = MainWindow::GetMainWindow()->GetLayerCollection("MRI");
  QList<Layer*> layers = col_mri->GetLayers();
  int nSel = 0;
  for ( int i = 0; i < layers.size(); i++ )
  {
    ui->comboBoxTemplate->addItem(layers[i]->GetName(), QVariant::fromValue((QObject*)layers[i]));
    if ( layers[i] == col_mri->GetActiveLayer() )
    {
      nSel = i;
    }
  }
  ui->comboBoxTemplate->setCurrentIndex(nSel);
  ui->lineEditName->setFocus();
}

DialogNewPointSet::~DialogNewPointSet()
{
  delete ui;
}

QString DialogNewPointSet::GetPointSetName()
{
  return ui->lineEditName->text().trimmed();
}

void DialogNewPointSet::SetPointSetName( const QString& name )
{
  ui->lineEditName->setText( name );
}

LayerMRI* DialogNewPointSet::GetTemplate()
{
  return qobject_cast<LayerMRI*>(
        ui->comboBoxTemplate->itemData(ui->comboBoxTemplate->currentIndex()).value<QObject*>());
}

void DialogNewPointSet::OnOK()
{
  if ( GetPointSetName().isEmpty())
  {
    QMessageBox::warning( this, "Error", "Point set name can not be empty." );
    return;
  }
  accept();
}

int DialogNewPointSet::GetType()
{
  return ( ui->radioButtonControlPoint->isChecked() ?
             LayerPropertyPointSet::ControlPoint :
             LayerPropertyPointSet::WayPoint );
}
