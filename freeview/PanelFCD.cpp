#include "PanelFCD.h"
#include "ui_PanelFCD.h"
#include "LayerFCD.h"
#include "LayerPropertyFCD.h"

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
}

void PanelFCD::DoUpdateWidgets()
{
  BlockAllSignals( true );
  /*
  for ( int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++ )
  {
    QTreeWidgetItem* item = ui->treeWidgetLayers->topLevelItem( i );
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole).value<QObject*>() );
    if ( layer )
    {
      item->setCheckState( 0, (layer->IsVisible() ? Qt::Checked : Qt::Unchecked) );
    }
  }
  */

  LayerFCD* layer = GetCurrentLayer<LayerFCD*>();

  /*
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  */

  if ( layer )
  {
    ui->sliderSigma->setValue((int)(layer->GetProperty()->GetSigma()));
    ui->sliderThreshold->setValue((int)(layer->GetProperty()->GetThicknessThreshold()));
    ChangeLineEditNumber(ui->lineEditSigma, layer->GetProperty()->GetSigma());
    ChangeLineEditNumber(ui->lineEditThreshold, layer->GetProperty()->GetThicknessThreshold());
  }

  BlockAllSignals( false );
}

void PanelFCD::DoIdle()
{

}

void PanelFCD::OnSliderThresholdReleased()
{

}

void PanelFCD::OnSliderSigmaReleased()
{

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
