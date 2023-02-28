#include "PanelFCD.h"
#include "ui_PanelFCD.h"
#include <cstddef>
#include "LayerFCD.h"
#include "LayerPropertyFCD.h"
#include "MainWindow.h"
#include <QDebug>
#include <QFileDialog>

PanelFCD::PanelFCD(QWidget *parent) :
  PanelLayer("FCD", parent),
  ui(new Ui::PanelFCD)
{
  ui->setupUi(this);
}

PanelFCD::~PanelFCD()
{
  delete ui;
}


void PanelFCD::ConnectLayer(Layer *layer_in)
{
  PanelLayer::ConnectLayer( layer_in );

  LayerFCD* layer = qobject_cast<LayerFCD*>(layer_in);
  if ( !layer )
  {
    return;
  }
  LayerPropertyFCD* p = layer->GetProperty();
  connect(p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect(layer, SIGNAL(LabelsChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
  connect(layer, SIGNAL(StatusChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
}

void PanelFCD::DoUpdateWidgets()
{
  BlockAllSignals( true );

  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();

  if ( layer )
  {
    ui->sliderSigma->setValue((int)(layer->GetProperty()->GetSigma()));
    ui->sliderThreshold->setValue((int)(layer->GetProperty()->GetThicknessThreshold()*10));
    ui->sliderMinArea->setValue(layer->GetProperty()->GetMinArea());
    ui->sliderOpacity->setValue((int)(layer->GetProperty()->GetOpacity()*100));
    ui->spinBoxOpacity->setValue(layer->GetProperty()->GetOpacity());

    ChangeLineEditNumber(ui->lineEditSigma, layer->GetProperty()->GetSigma());
    ChangeLineEditNumber(ui->lineEditThreshold, layer->GetProperty()->GetThicknessThreshold());
    ChangeLineEditNumber(ui->lineEditMinArea, layer->GetProperty()->GetMinArea());

    ui->sliderSigma->setEnabled(!layer->IsBusy());
    ui->sliderThreshold->setEnabled(!layer->IsBusy());
    ui->sliderMinArea->setEnabled(!layer->IsBusy());
    ui->lineEditSigma->setEnabled(!layer->IsBusy());
    ui->lineEditThreshold->setEnabled(!layer->IsBusy());
    ui->lineEditMinArea->setEnabled(!layer->IsBusy());
    ui->pushButtonRecompute->setEnabled(!layer->IsBusy());

    if (!layer->IsBusy())
      UpdateLabelList(layer);
  }

  BlockAllSignals( false );
}

void PanelFCD::UpdateLabelList(LayerFCD *layer_in)
{
  LayerFCD* layer = layer_in;
  if (!layer)
    layer = GetCurrentLayer<LayerFCD*>();

  ui->treeWidgetLabels->clear();
  if (layer && layer->GetFCDData())
  {
    FCD_DATA* fcd = layer->GetFCDData();
    QList<bool> flags = layer->GetLabelVisibility();
    for (int i = 0; i < fcd->nlabels; i++)
    {
      //   if (fcd->labels[i]->n_points > 0)
      {
        QTreeWidgetItem* item = new QTreeWidgetItem( ui->treeWidgetLabels );
        item->setText(0, QString("%1   %2mm").arg(QString(fcd->label_names[i]))
                      .arg(fcd->labels[i]->avg_stat, 0, 'g', 4));
        item->setCheckState(0, flags[i]?Qt::Checked:Qt::Unchecked);
      }
    }
  }
}

void PanelFCD::DoIdle()
{

}

void PanelFCD::OnSliderOpacityChanged(int)
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
    layer->GetProperty()->SetOpacity(ui->sliderOpacity->value()/100.0);
}

void PanelFCD::OnSliderThresholdReleased()
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
    layer->GetProperty()->SetThicknessThreshold(ui->sliderThreshold->value()/10.0);
}

void PanelFCD::OnSliderSigmaReleased()
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
    layer->GetProperty()->SetSigma(ui->sliderSigma->value());
}

void PanelFCD::OnSliderMinAreaReleased()
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
    layer->GetProperty()->SetMinArea(ui->sliderMinArea->value());
}

void PanelFCD::OnSliderThresholdChanged(int)
{
  ChangeLineEditNumber(ui->lineEditThreshold, ui->sliderThreshold->value()/10.0);
}

void PanelFCD::OnSliderSigmaChanged(int)
{
  ChangeLineEditNumber(ui->lineEditSigma, ui->sliderSigma->value());
}

void PanelFCD::OnSliderMinAreaChanged(int)
{
  ChangeLineEditNumber(ui->lineEditMinArea, ui->sliderMinArea->value());
}

void PanelFCD::OnTextSigmaReturned()
{
  bool ok;
  double val = ui->lineEditSigma->text().trimmed().toDouble(&ok);
  if (ok)
  {
    LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
    if ( layer )
      layer->GetProperty()->SetSigma(val);
  }
}

void PanelFCD::OnTextThresholdReturned()
{
  bool ok;
  double val = ui->lineEditThreshold->text().trimmed().toDouble(&ok);
  if (ok)
  {
    LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
    if ( layer )
      layer->GetProperty()->SetThicknessThreshold(val);
  }
}

void PanelFCD::OnTextMinAreaReturned()
{
  bool ok;
  int val = ui->lineEditMinArea->text().trimmed().toInt(&ok);
  if (ok)
  {
    LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
    if ( layer )
      layer->GetProperty()->SetMinArea(val);
  }
}

void PanelFCD::OnLabelSelectionChanged()
{
  int n = ui->treeWidgetLabels->indexOfTopLevelItem(ui->treeWidgetLabels->currentItem());
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if (layer && n >= 0)
  {
    double pos[3];
    layer->GetLabelCentroidPosition(n, pos);
    MainWindow::GetMainWindow()->SetSlicePosition(pos);
  }
}

void PanelFCD::OnLabelItemChanged(QTreeWidgetItem *item)
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
  {
    int n = ui->treeWidgetLabels->indexOfTopLevelItem(item);
    layer->SetLabelVisible(n, item->checkState(0) == Qt::Checked);
  }
}

void PanelFCD::OnButtonRecompute()
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
  {
    layer->GetProperty()->blockSignals(true);
    double val;
    int nval;
    bool ok;
    val = ui->lineEditThreshold->text().trimmed().toDouble(&ok);
    if (ok && val != layer->GetProperty()->GetThicknessThreshold())
      layer->GetProperty()->SetThicknessThreshold(val);

    val = ui->lineEditSigma->text().trimmed().toDouble(&ok);
    if (ok && val != layer->GetProperty()->GetSigma())
      layer->GetProperty()->SetSigma(val);

    nval = ui->lineEditMinArea->text().trimmed().toInt(&ok);
    if (ok && nval != layer->GetProperty()->GetMinArea())
      layer->GetProperty()->SetMinArea(nval);

    layer->GetProperty()->blockSignals(false);
    layer->Recompute();
  }
}

void PanelFCD::OnButtonGotoContralateral()
{
  MainWindow::GetMainWindow()->GoToContralateralPoint();
}

void PanelFCD::OnButtonSaveFCDLabels()
{
  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();
  if ( layer )
  {
    QString dir = QFileDialog::getExistingDirectory(
          this,
          tr("Select a directory to save FCD labels to"),
          QDir::currentPath() );
    if(!dir.isNull())
    {
      layer->SaveFCDLabels(dir);
    }
  }
}
