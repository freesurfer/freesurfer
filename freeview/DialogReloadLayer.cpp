#include "DialogReloadLayer.h"
#include "ui_DialogReloadLayer.h"
#include "Layer.h"

DialogReloadLayer::DialogReloadLayer(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogReloadLayer)
{
  ui->setupUi(this);
  ui->checkBoxCloseFirst->setChecked(true);
}

DialogReloadLayer::~DialogReloadLayer()
{
  delete ui;
}

int DialogReloadLayer::Execute(const QList<Layer*>& layers)
{
  QString type = layers.first()->GetPrimaryType();
  if (type == "MRI")
    type = "Volume";
  ui->checkBoxCloseFirst->setText(tr("Close selected %1(s) first before reloading them").arg(type));
  setWindowTitle(tr("Reload %1(s)").arg(type));
  QStringList filenames;
  foreach (Layer* layer, layers)
    filenames << layer->GetFileName();
  ui->labelMessage->setText(tr("Reload following %1(s):\n\n%2\n").arg(type).arg(filenames.join("\n")));
  return exec();
}

bool DialogReloadLayer::GetCloseLayerFirst()
{
  return ui->checkBoxCloseFirst->isChecked();
}
