#include "DialogReloadLayer.h"
#include "ui_DialogReloadLayer.h"



DialogReloadLayer::DialogReloadLayer(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogReloadLayer)
{
    ui->setupUi(this);
}

DialogReloadLayer::~DialogReloadLayer()
{
    delete ui;
}

int DialogReloadLayer::Execute(const QString &layer_name, const QString &layer_type, const QString &fn)
{
  setWindowTitle(tr("Reload %1").arg(layer_type));
  ui->labelMessage->setText(tr("Reload %1 %2 from:\n\n%3\n").arg(layer_type).arg(layer_name).arg(fn));
  return exec();
}

bool DialogReloadLayer::GetCloseLayerFirst()
{
  return ui->checkBoxCloseFirst->isChecked();
}
