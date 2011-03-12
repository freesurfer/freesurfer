#include "DialogNewVolume.h"
#include "ui_DialogNewVolume.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include <QMessageBox>

DialogNewVolume::DialogNewVolume(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogNewVolume)
{
    ui->setupUi(this);

    LayerCollection* col_mri = MainWindow::GetMainWindow()->GetLayerCollection("MRI");
    QList<Layer*> layers = col_mri->GetLayers();
    int nSel = 0;
    for ( int i = 0; i < layers.size(); i++ )
    {
        ui->comboBoxTemplate->addItem(layers[i]->GetName(), QVariant::fromValue((QObject*)layers[i]));
        if ( layers[i] == col_mri->GetActiveLayer() )
            nSel = i;
    }
    ui->comboBoxTemplate->setCurrentIndex(nSel);
    ui->lineEditName->setFocus();
}

DialogNewVolume::~DialogNewVolume()
{
    delete ui;
}

QString DialogNewVolume::GetVolumeName()
{
  return ui->lineEditName->text().trimmed();
}

void DialogNewVolume::SetVolumeName( const QString& name )
{
  ui->lineEditName->setText( name );
}

bool DialogNewVolume::GetCopyVoxel()
{
  return ui->checkBoxCopyData->isChecked();
}

void DialogNewVolume::SetCopyVoxel( bool bVoxel )
{
  ui->checkBoxCopyData->setChecked( bVoxel );
}

LayerMRI* DialogNewVolume::GetTemplate()
{
  return qobject_cast<LayerMRI*>(
          ui->comboBoxTemplate->itemData(ui->comboBoxTemplate->currentIndex()).value<QObject*>());
}

int DialogNewVolume::GetDataType()
{
  if ( ui->comboBoxDataType->currentIndex() == ui->comboBoxDataType->count()-1 )
    return GetTemplate()->GetDataType();
  else
    return ui->comboBoxDataType->currentIndex();
}

void DialogNewVolume::OnOK()
{
    if ( GetVolumeName().isEmpty())
    {
        QMessageBox::warning( this, "Error", "Volume name can not be empty." );
        return;
    }
    accept();
}

void DialogNewVolume::OnToggleCopyVoxelData(bool bCopy)
{
    if (bCopy)
        ui->comboBoxDataType->setCurrentIndex(ui->comboBoxDataType->count()-1);
    ui->comboBoxDataType->setDisabled(bCopy);
}
