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
#include "PanelTrack.h"
#include "ui_PanelTrack.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "LayerTrack.h"
#include "MyUtils.h"
#include "LayerCollection.h"
#include "LayerPropertyTrack.h"
#include <QFileInfo>

PanelTrack::PanelTrack(QWidget *parent) :
  PanelLayer("Tract", parent),
  ui(new Ui::PanelTrack)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (mainwnd)
  {
    ui->toolbar->addAction(mainwnd->ui->actionLoadTrack);
    ui->toolbar->addAction(mainwnd->ui->actionCloseTrack);
  }

  m_widgetlistDirectionalColor << ui->labelDirectionScheme
                               << ui->comboBoxDirectionScheme
                               << ui->labelDirectionMapping
                               << ui->comboBoxDirectionMapping;
  m_widgetlistSolidColor << ui->labelSolidColor
                         << ui->colorPickerSolidColor;
  m_widgetlistScalarColor << ui->labelScalar << ui->labelScalarColor << ui->labelScalarMin << ui->labelScalarMax
                          << ui->comboBoxScalar << ui->comboBoxScalarColor << ui->lineEditScalarMin << ui->lineEditScalarMax
                          << ui->sliderScalarMin << ui->sliderScalarMax;
  connect(ui->pushButtonShowClusterMap, SIGNAL(clicked()), mainwnd, SLOT(ShowTractClusterMap()));
  connect(ui->sliderOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderOpacity(int)));
  connect(ui->lineEditOpacity, SIGNAL(textChanged(QString)), SLOT(OnLineEditOpacity(QString)));
  connect(ui->sliderScalarMin, SIGNAL(valueChanged(int)), SLOT(OnSliderScalarThreshold(int)));
  connect(ui->sliderScalarMax, SIGNAL(valueChanged(int)), SLOT(OnSliderScalarThreshold(int)));
  connect(ui->lineEditScalarMin, SIGNAL(textChanged(QString)), SLOT(OnLineEditScalarThreshold(QString)));
  connect(ui->lineEditScalarMax, SIGNAL(textChanged(QString)), SLOT(OnLineEditScalarThreshold(QString)));
}

PanelTrack::~PanelTrack()
{
  delete ui;
}

void PanelTrack::ConnectLayer(Layer *layer_in)
{
  PanelLayer::ConnectLayer( layer_in );

  LayerTrack* layer = qobject_cast<LayerTrack*>(layer_in);
  if ( !layer )
  {
    return;
  }
  LayerPropertyTrack* p = layer->GetProperty();
  connect(p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect(ui->comboBoxColorCode, SIGNAL(currentIndexChanged(int)), p, SLOT(SetColorCode(int)) );
  connect(ui->comboBoxDirectionScheme, SIGNAL(currentIndexChanged(int)), p, SLOT(SetDirectionScheme(int)));
  connect(ui->comboBoxDirectionMapping, SIGNAL(currentIndexChanged(int)), p, SLOT(SetDirectionMapping(int)));
  connect(ui->colorPickerSolidColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetSolidColor(QColor)));
  connect(ui->comboBoxRenderRep, SIGNAL(currentIndexChanged(int)), p, SLOT(SetRenderRep(int)));
  connect(ui->comboBoxScalar, SIGNAL(currentIndexChanged(int)), p, SLOT(SetScalarIndex(int)));
  connect(ui->comboBoxScalarColor, SIGNAL(currentIndexChanged(int)), p, SLOT(SetScalarColorMap(int)));
}

void PanelTrack::DoUpdateWidgets()
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

  LayerTrack* layer = GetCurrentLayer<LayerTrack*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }

  ui->lineEditFileName->clear();
  ShowWidgets(m_widgetlistDirectionalColor, layer && layer->GetProperty()->GetColorCode() == LayerPropertyTrack::Directional);
  ShowWidgets(m_widgetlistSolidColor, layer && layer->GetProperty()->GetColorCode() == LayerPropertyTrack::SolidColor);
  ShowWidgets(m_widgetlistScalarColor, layer && layer->GetProperty()->GetColorCode() == LayerPropertyTrack::Scalar);
  if ( layer )
  {
    QString fn = layer->GetFileName();
    if (layer->IsCluster())
      fn = QFileInfo(fn).absolutePath() + "/*.trk";
    ui->lineEditFileName->setText(fn);
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
    ui->comboBoxColorCode->setCurrentIndex(layer->GetProperty()->GetColorCode());
    ui->comboBoxColorCode->setItemData(LayerPropertyTrack::Scalar, layer->GetScalarNames().isEmpty()?0:33, Qt::UserRole-1);
    ui->comboBoxColorCode->setItemData(LayerPropertyTrack::EmbeddedColor, layer->HasEmbeddedColor()?33:0, Qt::UserRole-1);
    ui->comboBoxDirectionMapping->setCurrentIndex(layer->GetProperty()->GetDirectionMapping());
    ui->comboBoxDirectionScheme->setCurrentIndex(layer->GetProperty()->GetDirectionScheme());
    ui->colorPickerSolidColor->setCurrentColor(layer->GetProperty()->GetSolidColor());
    ui->comboBoxRenderRep->setCurrentIndex(layer->GetProperty()->GetRenderRep());
    ChangeLineEditNumber(ui->lineEditOpacity, layer->GetProperty()->GetOpacity());
    ui->sliderOpacity->setValue(layer->GetProperty()->GetOpacity()*100);

    ui->comboBoxScalar->clear();
    QStringList scalarNames = layer->GetScalarNames();
    foreach (QString name, scalarNames)
      ui->comboBoxScalar->addItem(name);
    ui->comboBoxScalar->setCurrentIndex(layer->GetProperty()->GetScalarIndex());

    if (!scalarNames.isEmpty())
    {
      ui->comboBoxScalar->setCurrentIndex(layer->GetProperty()->GetScalarIndex());
      ui->comboBoxScalarColor->setCurrentIndex(layer->GetProperty()->GetScalarColorMap());
      double range[2], th[2];
      layer->GetProperty()->GetScalarThreshold(th);
      layer->GetScalarRange(range);
      if (range[1] <= range[0])
        range[1] = range[0]+1;
      ChangeLineEditNumber(ui->lineEditScalarMin, th[0]);
      ChangeLineEditNumber(ui->lineEditScalarMax, th[1]);
      ui->sliderScalarMin->setValue((th[0]-range[0])/(range[1]-range[0])*100);
      ui->sliderScalarMax->setValue((th[1]-range[0])/(range[1]-range[0])*100);
    }
  }

  ui->labelFileName->setEnabled( layer );
  ui->lineEditFileName->setEnabled( layer );
  ui->pushButtonShowClusterMap->setVisible(layer && layer->IsCluster());

  BlockAllSignals( false );
}

void PanelTrack::DoIdle()
{

}

void PanelTrack::OnSliderOpacity(int val)
{
  LayerTrack* layer = GetCurrentLayer<LayerTrack*>();
  if (layer)
  {
    layer->GetProperty()->SetOpacity(val/100.0);
  }
  ChangeLineEditNumber(ui->lineEditOpacity, val/100.0);
}

void PanelTrack::OnLineEditOpacity(const QString & text)
{
  LayerTrack* layer = GetCurrentLayer<LayerTrack*>();
  bool bOK;
  double val = text.toDouble(&bOK);
  if (layer && bOK)
  {
    layer->GetProperty()->SetOpacity(val);
    ui->sliderOpacity->blockSignals(true);
    ui->sliderOpacity->setValue(val*100);
    ui->sliderOpacity->blockSignals(false);
  }
}


void PanelTrack::OnSliderScalarThreshold(int val)
{
  LayerTrack* layer = GetCurrentLayer<LayerTrack*>();
  if (layer)
  {
    double range[2], th[2];
    layer->GetScalarRange(range);
    layer->GetProperty()->GetScalarThreshold(th);
    int n = (sender() == ui->sliderScalarMax)?1:0;
    th[n] = (range[1]-range[0])*val/100+range[0];
    layer->GetProperty()->SetScalarThreshold(th[0], th[1]);
    ChangeLineEditNumber(n == 0?ui->lineEditScalarMin:ui->lineEditScalarMax, th[n]);
  }
}

void PanelTrack::OnLineEditScalarThreshold(const QString & text)
{
  LayerTrack* layer = GetCurrentLayer<LayerTrack*>();
  bool bOK;
  double val = text.toDouble(&bOK);
  if (layer && bOK)
  {
    double range[2], th[2];
    if (range[1] <= range[0])
      range[1] = range[0]+1;
    layer->GetScalarRange(range);
    layer->GetProperty()->GetScalarThreshold(th);
    int n = (sender() == ui->lineEditScalarMax)?1:0;
    th[n] = val;
    layer->GetProperty()->SetScalarThreshold(th[0], th[1]);
    QSlider* slider = (n == 0)?ui->sliderScalarMin:ui->sliderScalarMax;
    slider->blockSignals(true);
    slider->setValue((val-range[0])/(range[1]-range[0])*100);
    slider->blockSignals(false);
  }
}
