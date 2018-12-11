#include "DialogSurfaceLabelOperations.h"
#include "ui_DialogSurfaceLabelOperations.h"

DialogSurfaceLabelOperations::DialogSurfaceLabelOperations(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSurfaceLabelOperations)
{
  ui->setupUi(this);
  connect(ui->pushButtonClose, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonOpen, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonDilate, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
  connect(ui->pushButtonErode, SIGNAL(clicked(bool)), SLOT(OnButtonClicked()));
}

DialogSurfaceLabelOperations::~DialogSurfaceLabelOperations()
{
  delete ui;
}

void DialogSurfaceLabelOperations::OnButtonClicked()
{
  QVariantMap op;
  if (sender() == ui->pushButtonDilate)
  {
    op["operation"] = "dilate";
    op["times"] = ui->spinBoxDilateTimes->value();
  }
  else if (sender() == ui->pushButtonErode)
  {
    op["operation"] = "erode";
    op["times"] = ui->spinBoxErodeTimes->value();
  }
  else if (sender() == ui->pushButtonOpen)
  {
    op["operation"] = "open";
    op["times"] = ui->spinBoxOpenTimes->value();
  }
  else if (sender() == ui->pushButtonClose)
  {
    op["operation"] = "close";
    op["times"] = ui->spinBoxCloseTimes->value();
  }

  if (!op.isEmpty())
    emit OperationTriggered(op);
}
