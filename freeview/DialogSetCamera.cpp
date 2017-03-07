#include "DialogSetCamera.h"
#include "ui_DialogSetCamera.h"
#include "MainWindow.h"
#include "RenderView3D.h"

DialogSetCamera::DialogSetCamera(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSetCamera)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint |
                        Qt::WindowCloseButtonHint| Qt::CustomizeWindowHint );
  ui->labelMessage->hide();
}

DialogSetCamera::~DialogSetCamera()
{
  delete ui;
}

void DialogSetCamera::OnRefresh()
{
  MainWindow::GetMainWindow()->GetRenderView(3)->Reset();
  QStringList cmd;
  cmd << "camera";
  cmd << "Azimuth" << QString::number(ui->doubleSpinBoxAzimuth->value());
  cmd << "Zoom" << QString::number(ui->doubleSpinBoxZoom->value());
  cmd << "Elevation" << QString::number(ui->doubleSpinBoxElevation->value());
  cmd << "Roll" << QString::number(ui->doubleSpinBoxRoll->value());
  MainWindow::GetMainWindow()->CommandSetCamera(cmd);
}

void DialogSetCamera::OnReset()
{
  QList<QWidget*> widgets = findChildren<QWidget*>();
  foreach(QWidget* w, widgets)
    w->blockSignals(true);

  ui->doubleSpinBoxAzimuth->setValue(0);
  ui->doubleSpinBoxElevation->setValue(0);
  ui->doubleSpinBoxRoll->setValue(0);
  ui->doubleSpinBoxZoom->setValue(1);

  foreach(QWidget* w, widgets)
    w->blockSignals(false);

  OnRefresh();
}

void DialogSetCamera::OnSettingChanged()
{
  OnRefresh();
}
