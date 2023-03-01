#include "DialogThresholdVolume.h"
#include "ui_DialogThresholdVolume.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include <QMessageBox>

DialogThresholdVolume::DialogThresholdVolume(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogThresholdVolume)
{
  ui->setupUi(this);
}

DialogThresholdVolume::~DialogThresholdVolume()
{
  delete ui;
}

void DialogThresholdVolume::OnButtonRun()
{
  if (ValidateInputs())
  {
    LayerMRI* src = qobject_cast<LayerMRI*>(ui->comboBoxSourceVolume->itemData(ui->comboBoxSourceVolume->currentIndex()).value<QObject*>());
    LayerMRI* target = qobject_cast<LayerMRI*>(ui->comboBoxTargetVolume->itemData(ui->comboBoxTargetVolume->currentIndex()).value<QObject*>());
    double th_high = ui->lineEditHighThreshold->text().toDouble();
    if (ui->lineEditHighThreshold->text().trimmed().toLower() == "inf")
      th_high = 1e20;
    target->Threshold(ui->lineEditTargetFrame->text().trimmed().toInt(), src,
                      ui->lineEditSourceFrame->text().trimmed().toInt(),
                      ui->lineEditLowThreshold->text().trimmed().toDouble(), th_high,
                      ui->checkBoxReplaceIn->isChecked(), ui->lineEditInValue->text().trimmed().toDouble(),
                      ui->checkBoxReplaceOut->isChecked(), ui->lineEditOutValue->text().trimmed().toDouble());
  }
}

void DialogThresholdVolume::OnButtonReset()
{
  LayerMRI* target = qobject_cast<LayerMRI*>(ui->comboBoxTargetVolume->itemData(ui->comboBoxTargetVolume->currentIndex()).value<QObject*>());
  if (target)
    target->RestoreFromBackup();
}

void DialogThresholdVolume::UpdateVolumes()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
  LayerMRI* src = NULL, *target = NULL;
  if ( ui->comboBoxSourceVolume->currentIndex() >= 0)
    src = qobject_cast<LayerMRI*>(ui->comboBoxSourceVolume->itemData(ui->comboBoxSourceVolume->currentIndex()).value<QObject*>());
  if ( ui->comboBoxTargetVolume->currentIndex() >= 0)
    target = qobject_cast<LayerMRI*>(ui->comboBoxTargetVolume->itemData(ui->comboBoxTargetVolume->currentIndex()).value<QObject*>());
  ui->comboBoxSourceVolume->clear();
  ui->comboBoxTargetVolume->clear();
  foreach (Layer* layer, layers)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
    ui->comboBoxSourceVolume->addItem(mri->GetName(), QVariant::fromValue<QObject*>(mri));
    if (src == mri)
      ui->comboBoxSourceVolume->setCurrentIndex(ui->comboBoxSourceVolume->count()-1);
    ui->comboBoxTargetVolume->addItem(mri->GetName(), QVariant::fromValue<QObject*>(mri));
    if (target == mri)
      ui->comboBoxTargetVolume->setCurrentIndex(ui->comboBoxTargetVolume->count()-1);
  }
}

bool DialogThresholdVolume::ValidateInputs()
{
  if (ui->comboBoxSourceVolume->currentIndex() < 0)
  {
    QMessageBox::warning(this, "Error", "Please select a source volume");
    return false;
  }
  LayerMRI* src = qobject_cast<LayerMRI*>(ui->comboBoxSourceVolume->itemData(ui->comboBoxSourceVolume->currentIndex()).value<QObject*>());
  LayerMRI* target = qobject_cast<LayerMRI*>(ui->comboBoxTargetVolume->itemData(ui->comboBoxTargetVolume->currentIndex()).value<QObject*>());
  bool ok;
  int source_frame = ui->lineEditSourceFrame->text().toInt(&ok);
  if (!ok || source_frame < 0)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid source frame number.");
    return false;
  }
  else if (source_frame >= src->GetNumberOfFrames())
  {
    QMessageBox::warning(this, "Error", "Source frame number is out of bound. Please enter a valid number.");
    return false;
  }
  int target_frame = ui->lineEditTargetFrame->text().toInt(&ok);
  if (!ok || target_frame < 0)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid target frame number.");
    return false;
  }
  else if (target_frame >= target->GetNumberOfFrames())
  {
    QMessageBox::warning(this, "Error", "Target frame number is out of bound. Please enter a valid number.");
    return false;
  }

  double th = ui->lineEditLowThreshold->text().toDouble(&ok);
  if (!ok)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid low threshold.");
    return false;
  }
  th = ui->lineEditHighThreshold->text().trimmed().toDouble(&ok);
  if (ui->lineEditHighThreshold->text().trimmed().toLower() != "inf" && !ok)
  {
    QMessageBox::warning(this, "Error", "Please enter a valid high threshold.");
    return false;
  }
  return true;
}
