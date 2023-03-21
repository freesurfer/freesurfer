#include "PanelConnectomeMatrix.h"
#include "ui_PanelConnectomeMatrix.h"
#include "LayerConnectomeMatrix.h"
#include "LayerPropertyConnectomeMatrix.h"
#include "MainWindow.h"

PanelConnectomeMatrix::PanelConnectomeMatrix(QWidget *parent) :
  PanelLayer("CMAT", parent),
  ui(new Ui::PanelConnectomeMatrix),
  m_bColorTableDirty(true)
{
  ui->setupUi(this);
}

PanelConnectomeMatrix::~PanelConnectomeMatrix()
{
  delete ui;
}

void PanelConnectomeMatrix::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerConnectomeMatrix* layer = qobject_cast<LayerConnectomeMatrix*>(layer_in);
  if ( !layer )
  {
    return;
  }
  m_bColorTableDirty = true;

  LayerPropertyConnectomeMatrix* p = layer->GetProperty();
  connect(p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
  connect(ui->doubleSpinBoxFromOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetFromLabelOpacity(double)));
  connect(ui->doubleSpinBoxToOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetToLabelOpacity(double)));
  connect( ui->colorPickerSpline, SIGNAL(colorChanged(QColor)), p, SLOT(SetSplineColor(QColor)));
  /*
  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->checkBoxSmooth, SIGNAL(stateChanged(int)), p, SLOT(SetTextureSmoothing(int)) );
  connect( ui->checkBoxShowContour, SIGNAL(clicked(bool)), p, SLOT(SetShowAsContour(bool)) );
  connect( ui->checkBoxShowLabelContour, SIGNAL(clicked(bool)), p, SLOT(SetShowAsLabelContour(bool)) );
  connect( ui->sliderFrame, SIGNAL(valueChanged(int)), layer, SLOT(SetActiveFrameOneBase(int)) );
  connect( ui->spinBoxFrame, SIGNAL(valueChanged(int)), layer, SLOT(SetActiveFrameOneBase(int)) );
  connect( ui->checkBoxDisplayVector, SIGNAL(toggled(bool)), p, SLOT(SetDisplayVector(bool)) );
  connect( ui->checkBoxDisplayTensor, SIGNAL(toggled(bool)), p, SLOT(SetDisplayTensor(bool)) );
  connect( ui->comboBoxRenderObject, SIGNAL(currentIndexChanged(int)), p, SLOT(SetVectorRepresentation(int)) );
  connect( ui->comboBoxInversion, SIGNAL(currentIndexChanged(int)), p, SLOT(SetVectorInversion(int)) );
  connect( ui->checkBoxProjectionMap, SIGNAL(toggled(bool)), p, SLOT(SetShowProjectionMap(bool)));
  if ( layer->IsTypeOf( "DTI" ) )
    connect( ui->comboBoxDirectionCode, SIGNAL(currentIndexChanged(int)),
             qobject_cast<LayerDTI*>(layer)->GetProperty(), SLOT(SetDirectionCode(int)) );
  connect( layer, SIGNAL(ActiveFrameChanged(int)), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(ActiveFrameChanged(int)), this, SLOT(OnActiveFrameChanged(int)));
  connect( layer, SIGNAL(FillValueChanged(double)), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(LabelStatsReady()), this, SLOT(UpdateWidgets()));
  connect( ui->checkBoxClearBackground, SIGNAL(toggled(bool)), p, SLOT(SetClearBackground(bool)) );
  connect( ui->checkBoxClearHigher, SIGNAL(toggled(bool)), p, SLOT(SetHeatScaleClearHigh(bool)) );
  connect( ui->checkBoxTruncate, SIGNAL(toggled(bool)), p, SLOT(SetHeatScaleTruncate(bool)) );
  connect( ui->checkBoxInvert, SIGNAL(toggled(bool)), p, SLOT(SetHeatScaleInvert(bool)) );
  connect( ui->checkBoxShowOutline, SIGNAL(toggled(bool)), p, SLOT(SetShowLabelOutline(bool)) );
  connect( ui->checkBoxContourExtractAll, SIGNAL(toggled(bool)), p, SLOT(SetContourExtractAllRegions(bool)) );
  connect( ui->checkBoxUseColorMap, SIGNAL(toggled(bool)), p, SLOT(SetContourUseImageColorMap(bool)) );
  connect( ui->checkBoxShowInfo, SIGNAL(toggled(bool)), p, SLOT(SetShowInfo(bool)) );
  connect( ui->colorPickerContour, SIGNAL(colorChanged(QColor)), p, SLOT(SetContourColor(QColor)));
  connect( ui->checkBoxUpsampleContour, SIGNAL(toggled(bool)), p, SLOT(SetContourUpsample(bool)));
  connect( ui->checkBoxRememberFrame, SIGNAL(toggled(bool)), p, SLOT(SetRememberFrameSettings(bool)));
  */
}

void PanelConnectomeMatrix::DoIdle()
{

}

void PanelConnectomeMatrix::DoUpdateWidgets()
{
  BlockAllSignals(true);

  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (layer)
  {
    ui->lineEditFileName->setText(layer->GetFileName());
    ui->lineEditFileName->setCursorPosition(ui->lineEditFileName->text().size());
    double dval = layer->GetProperty()->GetFromLabelOpacity();
    ui->sliderFromOpacity->setValue((int)(dval*100));
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxFromOpacity, dval );
    dval = layer->GetProperty()->GetToLabelOpacity();
    ui->sliderToOpacity->setValue((int)(dval*100));
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxToOpacity, dval );

    ChangeLineEditNumber(ui->lineEditSplineRadius, layer->GetProperty()->GetSplineRadius());
    ui->colorPickerSpline->setCurrentColor(layer->GetProperty()->GetSplineColor());
  }

  if (m_bColorTableDirty)
    PopulateColorTable();

  BlockAllSignals(false);
}

void PanelConnectomeMatrix::PopulateColorTable()
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (!layer)
    return;

  ui->treeWidgetFrom->clear();
  ui->treeWidgetTo->clear();
  QList<int> labels = layer->GetLabelList();
  foreach (int n, labels)
  {
    QString name = layer->GetLabelName(n);
    QColor color = layer->GetLabelColor(n);
    AddColorTableItem(n, name, color, ui->treeWidgetFrom);
    AddColorTableItem(n, name, color, ui->treeWidgetTo);
  }
  m_bColorTableDirty = false;
}

void PanelConnectomeMatrix::UpdateToLabelVisibility()
{
  QList<QTreeWidgetItem*> items = ui->treeWidgetFrom->selectedItems();
  QList<int> indices;
  foreach (QTreeWidgetItem* item, items)
    indices << ui->treeWidgetFrom->indexOfTopLevelItem(item);
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (indices.isEmpty() || !layer)
    return;

  for (int i = 0; i < ui->treeWidgetTo->topLevelItemCount(); i++)
  {
    bool bConnected = false;
    for (int j = 0; j < indices.size(); j++)
    {
      if (layer->HasConnection(indices[j], i))
      {
        bConnected = true;
        break;
      }
    }
    ui->treeWidgetTo->topLevelItem(i)->setHidden(!bConnected);
  }
}

void PanelConnectomeMatrix::AddColorTableItem(int value, const QString& name, const QColor& color,
                                              QTreeWidget *treeWidget)
{
  QTreeWidgetItem* item = new QTreeWidgetItem( treeWidget );
  item->setText(0, QString("%1 %2").arg(value).arg(name));
  QPixmap pix(13, 13);
  pix.fill( color );
  item->setIcon(0, QIcon(pix) );
  item->setData( 0, Qt::UserRole, color );
  item->setData(0, Qt::UserRole+1, value);
}

void PanelConnectomeMatrix::OnCurrentFromChanged()
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (!layer)
    return;

  QList<QTreeWidgetItem*> items = ui->treeWidgetFrom->selectedItems();
  QList<int> indices;
  foreach (QTreeWidgetItem* item, items)
    indices << ui->treeWidgetFrom->indexOfTopLevelItem(item);
  layer->SetFromLabelIndices(indices);

  UpdateToLabelVisibility();
}

void PanelConnectomeMatrix::OnCurrentToChanged()
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (!layer)
    return;

  QList<QTreeWidgetItem*> items = ui->treeWidgetTo->selectedItems();
  QList<int> indices;
  foreach (QTreeWidgetItem* item, items)
    indices << ui->treeWidgetTo->indexOfTopLevelItem(item);
  layer->SetToLabelIndices(indices);
}

void PanelConnectomeMatrix::OnCheckBoxToAll(bool bChecked)
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (!layer)
    return;

  layer->SetToAllLabels(bChecked);
}

void PanelConnectomeMatrix::OnSliderFromOpacity(int val)
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (layer)
  {
    layer->GetProperty()->SetFromLabelOpacity(val/100.0);
  }
}

void PanelConnectomeMatrix::OnSliderToOpacity(int val)
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  if (layer)
  {
    layer->GetProperty()->SetToLabelOpacity(val/100.0);
  }
}

void PanelConnectomeMatrix::OnLineEditSplineRadius(const QString &strg)
{
  LayerConnectomeMatrix* layer = GetCurrentLayer<LayerConnectomeMatrix*>();
  bool bOK;
  double val = strg.toDouble(&bOK);
  if (layer && bOK && val > 0)
  {
    layer->GetProperty()->SetSplineRadius(val);
  }
}
