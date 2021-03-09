/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "PanelROI.h"
#include "ui_PanelROI.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QToolBar>
#include "LayerROI.h"
#include "LayerPropertyROI.h"
#include "MyUtils.h"
#include "LayerSurface.h"
#include <QDebug>

PanelROI::PanelROI(QWidget *parent) :
  PanelLayer("ROI", parent),
  ui(new Ui::PanelROI)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (mainwnd)
  {
    ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionNewROI);
    ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionLoadROI);
    ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionCloseROI);
    ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionSaveROI);
    ui->toolbar->insertSeparator(ui->actionMoveLayerUp);

    m_listWidgetsHeatscale << ui->labelHeatscaleRange
                           << ui->lineEditHeatscaleMin << ui->lineEditHeatscaleMax;
    connect(mainwnd->GetLayerCollection("Surface"), SIGNAL(LayerAdded(Layer*)),
            this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
    connect(mainwnd->GetLayerCollection("Surface"), SIGNAL(LayerRemoved(Layer*)),
            this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
  }
}

PanelROI::~PanelROI()
{
  delete ui;
}

void PanelROI::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerROI* layer = qobject_cast<LayerROI*>(layer_in);
  if ( !layer )
  {
    return;
  }

  LayerPropertyROI* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->colorPickerColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetColor(QColor)) );
  connect( ui->comboBoxColorMap, SIGNAL(currentIndexChanged(int)), p, SLOT(SetColorCode(int)));
}

void PanelROI::DoIdle()
{
  // update action status
  BlockAllSignals( true );

  BlockAllSignals( false );
}

void PanelROI::OnSliderOpacity(int nVal)
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if ( layer )
  {
    layer->GetProperty()->SetOpacity( nVal / 100.0 );
  }
}

void PanelROI::OnEditThreshold(const QString& text)
{
  bool ok;
  double th = text.trimmed().toDouble(&ok);
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if ( ok && layer )
  {
    layer->GetProperty()->SetThreshold(th);
  }
}

void PanelROI::OnEditHeatscaleMin(const QString &text)
{
  bool ok;
  double val = text.trimmed().toDouble(&ok);
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if ( ok && layer )
  {
    layer->GetProperty()->SetHeatscaleMin(val);
  }
}

void PanelROI::OnEditHeatscaleMax(const QString &text)
{
  bool ok;
  double val = text.trimmed().toDouble(&ok);
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if ( ok && layer )
  {
    layer->GetProperty()->SetHeatscaleMax(val);
  }
}

void PanelROI::DoUpdateWidgets()
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

  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }

  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );

    double* rgb = layer->GetProperty()->GetColor();
    ui->colorPickerColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );

    ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );

    ChangeLineEditNumber(ui->lineEditThreshold, layer->GetProperty()->GetThreshold());

    ui->comboBoxColorMap->setCurrentIndex(layer->GetProperty()->GetColorCode());
    ShowWidgets(m_listWidgetsHeatscale, layer->GetProperty()->GetColorCode() == LayerPropertyROI::Heatscale);
    ui->colorPickerColor->setVisible(layer->GetProperty()->GetColorCode() == LayerPropertyROI::SolidColor);

    ChangeLineEditNumber(ui->lineEditHeatscaleMin, layer->GetProperty()->GetHeatscaleMin());
    ChangeLineEditNumber(ui->lineEditHeatscaleMax, layer->GetProperty()->GetHeatscaleMax());
  }

  ui->doubleSpinBoxOpacity->setEnabled( layer );
  ui->colorPickerColor->setEnabled( layer );
  ui->sliderOpacity->setEnabled( layer );
  ui->labelColor->setEnabled( layer );
  ui->labelFileName->setEnabled( layer );
  ui->lineEditFileName->setEnabled( layer );
  ui->labelOpacity->setEnabled( layer );
  ui->comboBoxMappedSurface->setEnabled(layer);

  QList<Layer*> surfs = MainWindow::GetMainWindow()->GetLayers("Surface");
  ui->comboBoxMappedSurface->clear();
  ui->comboBoxMappedSurface->addItem("None");
  int nIndex = 0;
  if (layer)
  {
    for (int i = 0; i < surfs.size(); i++)
    {
      ui->comboBoxMappedSurface->addItem(surfs[i]->GetName(), QVariant::fromValue((QObject*)surfs[i]));
      if (surfs[i] == layer->GetMappedSurface())
        nIndex = i+1;
    }
  }
  ui->comboBoxMappedSurface->setCurrentIndex(nIndex);
//  ui->widgetDilateErode->setVisible(nIndex > 0);
  ui->pushButtonResample->setVisible(nIndex > 0);

  BlockAllSignals( false );
}

void PanelROI::OnComboMappedSurface(int nIndex)
{
  LayerSurface* surf = qobject_cast<LayerSurface*>(ui->comboBoxMappedSurface->itemData(nIndex).value<QObject*>());
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
  {
    layer->SetMappedSurface(surf);
    UpdateWidgets();
  }
}

void PanelROI::OnButtonDilate()
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
    layer->Dilate(ui->spinBoxDilateTimes->value());
}

void PanelROI::OnButtonErode()
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
    layer->Erode(ui->spinBoxErodeTimes->value());
}

void PanelROI::OnButtonOpen()
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
    layer->Open(ui->spinBoxOpenTimes->value());
}

void PanelROI::OnButtonClose()
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
    layer->Close(ui->spinBoxCloseTimes->value());
}

void PanelROI::OnButtonResample()
{
  LayerROI* layer = GetCurrentLayer<LayerROI*>();
  if (layer)
    layer->Resample();
}
