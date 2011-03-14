#include "DialogGradientFilter.h"
#include "ui_DialogGradientFilter.h"

DialogGradientFilter::DialogGradientFilter(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogGradientFilter)
{
    ui->setupUi(this);
}

DialogGradientFilter::~DialogGradientFilter()
{
    delete ui;
}

void DialogGradientFilter::SetSmoothing(bool smooth)
{
    ui->checkBoxSmooth->setChecked(smooth);
}

bool DialogGradientFilter::GetSmoothing()
{
    return ui->checkBoxSmooth->isChecked();
}

void DialogGradientFilter::SetSD(double val)
{
    ui->doubleSpinBoxSD->setValue(val);
}

double DialogGradientFilter::GetSD()
{
    return ui->doubleSpinBoxSD->value();
}
