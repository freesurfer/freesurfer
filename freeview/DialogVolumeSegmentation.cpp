#include "DialogVolumeSegmentation.h"
#include "ui_DialogVolumeSegmentation.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include <QDebug>
#include <QMessageBox>

DialogVolumeSegmentation::DialogVolumeSegmentation(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogVolumeSegmentation)
{
  ui->setupUi(this);
}

DialogVolumeSegmentation::~DialogVolumeSegmentation()
{
  delete ui;
}

void DialogVolumeSegmentation::OnButtonRestore()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(ui->comboBoxVolume->itemData(ui->comboBoxVolume->currentIndex()).value<QObject*>());
  if (mri)
    mri->RestoreFromBackup();
}

void DialogVolumeSegmentation::OnButtonRun()
{
  if (!ValidateInput())
    return;

  LayerMRI* mri = qobject_cast<LayerMRI*>(ui->comboBoxVolume->itemData(ui->comboBoxVolume->currentIndex()).value<QObject*>());
  mri->Segment(ui->lineEditMinLabelIndex->text().trimmed().toInt(),
               ui->lineEditMaxLabelIndex->text().trimmed().toInt(),
               ui->lineEditMinNumberOfVoxels->text().trimmed().toInt());
}

bool DialogVolumeSegmentation::ValidateInput()
{
  if (ui->comboBoxVolume->currentIndex() < 0)
  {
    QMessageBox::warning(this, "Error", "Please select a volume.");
    return false;
  }
  bool ok;
  int nval = ui->lineEditMinLabelIndex->text().trimmed().toInt(&ok);
  if (!ok || nval < 0)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid minimum label index");
    return false;
  }
  nval = ui->lineEditMaxLabelIndex->text().trimmed().toInt(&ok);
  if (!ok || nval < 0)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid maximum label index");
    return false;
  }
  nval = ui->lineEditMinNumberOfVoxels->text().trimmed().toInt(&ok);
  if (!ok || nval < 0)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid minimum number of voxels");
    return false;
  }
  return true;
}

void DialogVolumeSegmentation::UpdateVolumes()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
  LayerMRI* src = NULL;
  if ( ui->comboBoxVolume->currentIndex() >= 0)
    src = qobject_cast<LayerMRI*>(ui->comboBoxVolume->itemData(ui->comboBoxVolume->currentIndex()).value<QObject*>());
  ui->comboBoxVolume->clear();
  foreach (Layer* layer, layers)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
    ui->comboBoxVolume->addItem(mri->GetName(), QVariant::fromValue<QObject*>(mri));
    if (src == mri)
      ui->comboBoxVolume->setCurrentIndex(ui->comboBoxVolume->count()-1);
  }
}
