#include "DialogCustomFill.h"
#include "ui_DialogCustomFill.h"
#include <QSettings>

DialogCustomFill::DialogCustomFill(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::DialogCustomFill)
{
  ui->setupUi(this);
  this->setWindowFlags(Qt::Tool);
  QSettings s;
  ui->checkBoxFillToPaths->setChecked(s.value("CustomFill/FillToPaths", true).toBool());
  ui->checkBoxFillToLabels->setChecked(s.value("CustomFill/FillToLabels").toBool());
  ui->checkBoxFillToThreashold->setChecked(s.value("CustomFill/FillToThreshold", true).toBool());
  ui->checkBoxFillToCurvature->setChecked(s.value("CustomFill/FillToCurvature").toBool());
  QVariant v = s.value("DialogCustomFill/Geometry");
  if (v.isValid())
  {
    this->restoreGeometry(v.toByteArray());
  }
}

DialogCustomFill::~DialogCustomFill()
{
  QSettings s;
  s.setValue("CustomFill/FillToPaths", ui->checkBoxFillToPaths->isChecked());
  s.setValue("CustomFill/FillToLabels", ui->checkBoxFillToLabels->isChecked());
  s.setValue("CustomFill/FillToThreshold", ui->checkBoxFillToThreashold->isChecked());
  s.setValue("CustomFill/FillToCurvature", ui->checkBoxFillToCurvature->isChecked());
  s.setValue("DialogCustomFill/Geometry", this->saveGeometry());

  delete ui;
}

void DialogCustomFill::OnButtonFill()
{
  QVariantMap map;
  map["DoNotMarkSurface"] = true;
  map["DoNotCrossPaths"] = ui->checkBoxFillToPaths->isChecked();
  map["DoNotCrossLabels"] = ui->checkBoxFillToLabels->isChecked();
  map["DoNotCrossThreshold"] = ui->checkBoxFillToThreashold->isChecked();
  map["DoNotFillUnlabeled"] = ui->checkBoxFillToUnlabeled->isChecked();
  map["FillToCurvature"] = ui->checkBoxFillToCurvature->isChecked();
  map["CreateLabel"] = ui->radioButtonCreateLabel->isChecked();
  map["AddToLabel"] = ui->radioButtonAddToLabel->isChecked();
  map["RemoveFromLabel"] = ui->radioButtonRemoveFromLabel->isChecked();
  map["UseAllPoints"] = ui->checkBoxUseAllPoints->isChecked();
  emit CustomFillTriggered(map);
}
