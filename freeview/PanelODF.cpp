#include "PanelODF.h"
#include "ui_PanelODF.h"
#include "LayerPropertyODF.h"
#include "LayerODF.h"

PanelODF::PanelODF(QWidget *parent) :
  PanelLayer("ODF", parent),
  ui(new Ui::PanelODF)
{
  ui->setupUi(this);
}

PanelODF::~PanelODF()
{
  delete ui;
}

void PanelODF::ConnectLayer(Layer *layer_in)
{
  PanelLayer::ConnectLayer( layer_in );

  LayerODF* layer = qobject_cast<LayerODF*>(layer_in);
  if ( !layer )
  {
    return;
  }
  LayerPropertyODF* p = layer->GetProperty();
  connect(p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect(layer, SIGNAL(LabelsChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
  connect(layer, SIGNAL(StatusChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
}

void PanelODF::DoUpdateWidgets()
{
  BlockAllSignals( true );

  LayerODF* layer = GetCurrentLayer<LayerODF*>();

  if ( layer )
  {
//    ui->sliderSigma->setValue((int)(layer->GetProperty()->GetSigma()));
//    ui->sliderThreshold->setValue((int)(layer->GetProperty()->GetThicknessThreshold()*10));
//    ui->sliderMinArea->setValue(layer->GetProperty()->GetMinArea());
//    ui->sliderOpacity->setValue((int)(layer->GetProperty()->GetOpacity()*100));
//    ui->spinBoxOpacity->setValue(layer->GetProperty()->GetOpacity());

//    ChangeLineEditNumber(ui->lineEditSigma, layer->GetProperty()->GetSigma());
//    ChangeLineEditNumber(ui->lineEditThreshold, layer->GetProperty()->GetThicknessThreshold());
//    ChangeLineEditNumber(ui->lineEditMinArea, layer->GetProperty()->GetMinArea());

//    ui->sliderSigma->setEnabled(!layer->IsBusy());
//    ui->sliderThreshold->setEnabled(!layer->IsBusy());
//    ui->sliderMinArea->setEnabled(!layer->IsBusy());
//    ui->lineEditSigma->setEnabled(!layer->IsBusy());
//    ui->lineEditThreshold->setEnabled(!layer->IsBusy());
//    ui->lineEditMinArea->setEnabled(!layer->IsBusy());
//    ui->pushButtonRecompute->setEnabled(!layer->IsBusy());
  }

  BlockAllSignals( false );
}

void PanelODF::DoIdle()
{

}
