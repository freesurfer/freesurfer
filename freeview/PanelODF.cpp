#include "PanelODF.h"
#include "ui_PanelODF.h"
#include "LayerPropertyODF.h"
#include "LayerODF.h"
#include "MainWindow.h"

PanelODF::PanelODF(QWidget *parent) :
  PanelLayer("ODF", parent),
  ui(new Ui::PanelODF)
{
  ui->setupUi(this);

  m_colorThresholdWidgets << ui->labelColorMin << ui->labelColorMax
                          << ui->lineEditColorMin << ui->lineEditColorMax
                          << ui->sliderColorMin << ui->sliderColorMax;
  connect(ui->comboBoxMask, SIGNAL(currentIndexChanged(int)), this, SLOT(OnComboMask(int)), Qt::QueuedConnection);
  connect(ui->lineEditScale, SIGNAL(textChanged(QString)), SLOT(OnLineEditScale(QString)), Qt::QueuedConnection);
  connect(ui->lineEditMaskThreshold, SIGNAL(returnPressed()), SLOT(OnLineEditMaskThresholdChanged()), Qt::QueuedConnection);
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
  connect(p, SIGNAL(OdfPropertyChanged()), this, SLOT(UpdateWidgets()), Qt::QueuedConnection );
  connect(p, SIGNAL(ColorCodeChanged()), this, SLOT(UpdateWidgets()), Qt::QueuedConnection );
  connect(ui->comboBoxInversion, SIGNAL(currentIndexChanged(int)), p, SLOT(SetOdfInversion(int)));
  connect(ui->spinBoxSkip, SIGNAL(valueChanged(int)), p, SLOT(SetOdfSkip(int)), Qt::QueuedConnection);
  connect(ui->comboBoxColorCode, SIGNAL(currentIndexChanged(int)), p, SLOT(SetOdfColorCode(int)));
  connect(ui->checkBoxShowIn2D, SIGNAL(toggled(bool)), p, SLOT(SetShowIn2DView(bool)), Qt::QueuedConnection);
}

void PanelODF::DoUpdateWidgets()
{
  BlockAllSignals( true );

  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->lineEditFileName->setText(layer->GetFileName());
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );

    ui->comboBoxMask->clear();
    ui->comboBoxMask->addItem("None");
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
    int n = 0;
    for (int i = 0; i < layers.size(); i++)
    {
      if (layer != layers[i])
      {
        ui->comboBoxMask->addItem( layers[i]->GetName(),  QVariant::fromValue((QObject*)layers[i]) );
        if (layer->GetOdfMask() == layers[i])
          n = ui->comboBoxMask->count()-1;
      }
    }
    ui->comboBoxMask->setCurrentIndex(n);

    ui->widgetMaskControl->setVisible(n > 0);
    double th[2];
    layer->GetOdfMaskThreshold(th);
    ChangeLineEditNumber(ui->lineEditMaskThreshold, th[0]);
    ChangeLineEditNumber(ui->lineEditScale, layer->GetProperty()->GetOdfScale());
    ui->spinBoxSkip->setValue(layer->GetProperty()->GetOdfSkip());
    ui->comboBoxInversion->setCurrentIndex( layer->GetProperty()->GetOdfInversion() );
    ui->comboBoxColorCode->setCurrentIndex( layer->GetProperty()->GetOdfColorCode());
    for (int i = 0; i < m_colorThresholdWidgets.size(); i++)
      m_colorThresholdWidgets[i]->setVisible(layer->GetProperty()->GetOdfColorCode() == 1);

    layer->GetProperty()->GetMagnitudeThreshold(th);
    double range_min = layer->GetProperty()->GetMinValue();
    double range_max = layer->GetProperty()->GetMaxValue();
    ChangeLineEditNumber( ui->lineEditColorMin, th[0] );
    ChangeLineEditNumber( ui->lineEditColorMax, th[1] );
    ui->sliderColorMin->setValue( (int)( ( th[0] - range_min ) / ( range_max - range_min ) * 100 ) );
    ui->sliderColorMax->setValue( (int)( ( th[1] - range_min ) / ( range_max - range_min ) * 100 ) );

    ui->checkBoxShowIn2D->setChecked(layer->GetProperty()->GetShowIn2DView());
  }

  BlockAllSignals( false );
}

void PanelODF::DoIdle()
{

}

void PanelODF::OnComboMask(int sel)
{
  LayerMRI* mask = qobject_cast<LayerMRI*>(ui->comboBoxMask->itemData(sel).value<QObject*>());
  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  if ( layer )
  {
    layer->SetOdfMask(mask);
  }
  ui->widgetMaskControl->setVisible(mask != NULL);
}

void PanelODF::OnLineEditMaskThresholdChanged()
{
  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  if (layer)
  {
    bool bOK;
    double dVal = ui->lineEditMaskThreshold->text().trimmed().toDouble( &bOK );
    if ( bOK )
    {
      double th[2];
      layer->GetOdfMaskThreshold(th);
      th[0] = dVal;
      layer->SetOdfMaskThreshold(th);
    }
  }
}

void PanelODF::OnLineEditScale( const QString& text )
{
  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  if (layer)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( bOK && dVal > 0)
    {
      layer->GetProperty()->SetOdfScale(dVal);
    }
  }
}

void PanelODF::OnLineEditColorThreshold(const QString &strg)
{
  bool bOk;
  double val = strg.toDouble(&bOk);
  if (!bOk)
    return;

  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  if (layer)
  {
    double th[2];
    layer->GetProperty()->GetMagnitudeThreshold(th);
    int n = (sender() == ui->lineEditColorMin ? 0:1);
    th[n] = val;
    layer->GetProperty()->SetMagnitudeThreshold(th);
    QSlider* slider = (n == 0?ui->sliderColorMin:ui->sliderColorMax);
    double range_min = layer->GetProperty()->GetMinValue();
    double range_max = layer->GetProperty()->GetMaxValue();
    slider->blockSignals(true);
    ui->sliderColorMin->setValue( (int)( ( th[n] - range_min ) / ( range_max - range_min ) * 100 ) );
    slider->blockSignals(false);
  }
}

void PanelODF::OnSliderColorThreshold(int nVal)
{
  LayerODF* layer = GetCurrentLayer<LayerODF*>();
  if (layer)
  {
    double th[2];
    layer->GetProperty()->GetMagnitudeThreshold(th);
    int n = (sender() == ui->sliderColorMin ? 0:1);
    double range_min = layer->GetProperty()->GetMinValue();
    double range_max = layer->GetProperty()->GetMaxValue();
    th[n] = range_min + (range_max-range_min)*nVal/100;
    layer->GetProperty()->SetMagnitudeThreshold(th);
    ChangeLineEditNumber(n == 0?ui->lineEditColorMin:ui->lineEditColorMax, th[n]);
  }
}
