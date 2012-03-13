#include "DialogLabelStats.h"
#include "ui_DialogLabelStats.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"

DialogLabelStats::DialogLabelStats(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DialogLabelStats)
{
    ui->setupUi(this);
    setWindowFlags(Qt::Dialog);
}

DialogLabelStats::~DialogLabelStats()
{
    delete ui;
}

void DialogLabelStats::OnSlicePositionChanged()
{
  if (!isVisible())
    return;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainwnd->GetLayers("MRI");
  LayerMRI* label = NULL, *mri = NULL;
  foreach (Layer* layer, layers)
  {
    LayerMRI* mri_layer = (LayerMRI*)layer;
    if (!label && mri_layer->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
      label = mri_layer;
    if (!mri &&  mri_layer->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT)
      mri = mri_layer;
    if (label && mri)
      break;
  }

  if (label)
  {
    float fLabel, fArea = 0;
    double mean, sd;
    int nCount = 0;
    label->GetCurrentLabelStats( mainwnd->GetMainViewId(), &fLabel, &nCount, &fArea, mri, &mean, &sd );
    ui->labelLabel->setText(QString::number((int)fLabel));
    ui->labelCount->setText(QString::number(nCount));
    ui->labelArea->setText(QString("%3 mm2").arg(fArea));
    if (mri)
      ui->labelMean->setText(QString("%1 +/- %2").arg(mean).arg(sd));
    else
      ui->labelMean->clear();
  }
}

void DialogLabelStats::showEvent(QShowEvent *e)
{
  OnSlicePositionChanged();
  QWidget::showEvent(e);
}
