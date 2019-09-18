#include "DialogLabelStats.h"
#include "ui_DialogLabelStats.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerROI.h"

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

void DialogLabelStats::UpdateStats()
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
  // only use current layer
  mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (mri->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
  {
    mri = NULL;
  }
  LayerROI* roi = (LayerROI*)mainwnd->GetActiveLayer("ROI");

  if (mri)
  {
    if (mri->GetNumberOfFrames() > 1)
      ui->labelVolume->setText(QString("%1 (frame %2)").arg(mri->GetName()).arg(mri->GetActiveFrame()));
    else
      ui->labelVolume->setText(mri->GetName());
  }
  else
    ui->labelVolume->clear();

  float fLabel, fArea = 0;
  double mean, sd;
  int nCount = 0;
  if (label && label->IsVisible())
  {
    label->GetCurrentLabelStats( mainwnd->GetMainViewId(), &fLabel, &nCount, &fArea, mri, &mean, &sd );
    ui->labelLabel->setText(QString::number((int)fLabel));
    ui->labelCount->setText(QString::number(nCount));
    ui->labelArea->setText(QString("%1 mm2").arg(fArea));
    if (mri)
    {
      ui->labelMean->setText(QString("%1 +/- %2").arg(mean).arg(sd));
      if (sd > 0)
        ui->labelSNR->setText(QString("%1").arg(mean/sd));
      else
        ui->labelSNR->setText("N/A");
    }
    else
    {
      ui->labelMean->clear();
      ui->labelSNR->clear();
    }
  }
  else if (roi) // update ROI
  {
    roi->GetStats( mainwnd->GetMainViewId(), &nCount, &fArea, mri, &mean, &sd );
    ui->labelLabel->setText(roi->GetName());
    ui->labelCount->setText(QString::number(nCount));
    ui->labelArea->setText("Volume:");
    ui->labelArea->setText(QString("%1 mm3").arg(fArea));
    if (mri)
    {
      ui->labelMean->setText(QString("%1 +/- %2").arg(mean).arg(sd));
      if (sd > 0)
        ui->labelSNR->setText(QString("%1").arg(mean/sd));
      else
        ui->labelSNR->setText("N/A");
    }
    else
    {
      ui->labelMean->clear();
      ui->labelSNR->clear();
    }
  }
}

void DialogLabelStats::showEvent(QShowEvent *e)
{
  UpdateStats();
  QWidget::showEvent(e);
}
