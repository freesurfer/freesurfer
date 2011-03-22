/**
 * @file  PanelPointSet.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/22 21:21:26 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "PanelPointSet.h"
#include "ui_PanelPointSet.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QToolBar>
#include "LayerPointSet.h"
#include "LayerPropertyPointSet.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <QFileDialog>
#include "MyUtils.h"

PanelPointSet::PanelPointSet(QWidget *parent) :
  PanelLayer(parent),
  ui(new Ui::PanelPointSet)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (mainwnd)
  {
    ui->toolbar->addAction(mainwnd->ui->actionNewPointSet);
    ui->toolbar->addAction(mainwnd->ui->actionLoadPointSet);
    ui->toolbar->addAction(mainwnd->ui->actionSavePointSet);
    ui->toolbar->addAction(mainwnd->ui->actionClosePointSet);
  }

  m_widgetlistSolidColor << ui->colorpickerSplineColor;

  m_widgetlistHeatScale << ui->comboBoxScalarMap
                        << ui->labelScalarMap
                        << ui->sliderMax
                        << ui->sliderMid
                        << ui->sliderMin
                        << ui->sliderOffset
                        << ui->lineEditMax
                        << ui->lineEditMid
                        << ui->lineEditMin
                        << ui->lineEditOffset
                        << ui->labelMax
                        << ui->labelMid
                        << ui->labelMin
                        << ui->labelOffset;

  m_widgetlistSpline << ui->labelSplineColor
                     << ui->comboBoxSplineColor
                     << ui->lineEditSplineRadius
                     << ui->labelSplineRadius;

  LayerCollection* lc = mainwnd->GetLayerCollection("PointSet");
  PanelLayer::InitializeLayerList( ui->treeWidgetLayers, lc );
}

PanelPointSet::~PanelPointSet()
{
  delete ui;
}


void PanelPointSet::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerPointSet* layer = qobject_cast<LayerPointSet*>(layer_in);
  if ( !layer )
  {
    return;
  }

  LayerPropertyPointSet* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->checkBoxShowSpline, SIGNAL(toggled(bool)), p, SLOT(SetShowSpline(bool)) );
  connect( ui->checkBoxSnapToCenter, SIGNAL(toggled(bool)), p, SLOT(SetSnapToVoxelCenter(bool)));
  connect( ui->colorpickerPointColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetColor(QColor)));
  connect( ui->colorpickerSplineColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetSplineColor(QColor)));
  connect( ui->comboBoxSplineColor, SIGNAL(currentIndexChanged(int)), p, SLOT(SetColorMap(int)));
}

void PanelPointSet::DoIdle()
{
  // update action status
  BlockAllSignals( true );


  BlockAllSignals( false );
}

void PanelPointSet::DoUpdateWidgets()
{
  BlockAllSignals( true );
  for ( int i = 0; i < ui->treeWidgetLayers->topLevelItemCount(); i++ )
  {
    QTreeWidgetItem* item = ui->treeWidgetLayers->topLevelItem( i );
    Layer* layer = qobject_cast<Layer*>( item->data(0, Qt::UserRole).value<QObject*>() );
    if ( layer )
    {
      item->setCheckState( 0, (layer->IsVisible() ? Qt::Checked : Qt::Unchecked) );
    }
  }

  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  int nColorMap = 0;
  bool bShowSpline = false;
  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );
    double* rgb = layer->GetProperty()->GetColor();
    ui->colorpickerPointColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetSplineColor();
    ui->colorpickerSplineColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
    ChangeLineEditNumber( ui->lineEditRadius, layer->GetProperty()->GetRadius() );
    ChangeLineEditNumber( ui->lineEditSplineRadius, layer->GetProperty()->GetSplineRadius() );

    nColorMap = layer->GetProperty()->GetColorMap();
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    ui->sliderMin->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMin() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderMid->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMid() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderMax->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMax() - fMin ) / ( fMax - fMin ) * 100 ) );
    ui->sliderOffset->setValue( (int)( ( layer->GetProperty()->GetHeatScaleOffset() + fMax ) / ( fMax + fMax ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditMid, layer->GetProperty()->GetHeatScaleMin() );
    ChangeLineEditNumber( ui->lineEditMin, layer->GetProperty()->GetHeatScaleMid() );
    ChangeLineEditNumber( ui->lineEditMax, layer->GetProperty()->GetHeatScaleMax() );
    ChangeLineEditNumber( ui->lineEditOffset, layer->GetProperty()->GetHeatScaleOffset() );

    ui->comboBoxSplineColor->setCurrentIndex( nColorMap );

    ui->comboBoxScalarMap->clear();;
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetLayers();
    int nSel = -1;
    for ( int i = 0; i < layers.size(); i++ )
    {
      ui->comboBoxScalarMap->addItem( layers[i]->GetName(), QVariant::fromValue((QObject*)layers[i]) );
      if ( layer->GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarLayer &&
           layer->GetProperty()->GetScalarLayer() == layers[i] )
      {
        nSel = i;
      }
    }
    std::vector<ScalarValues> svs = layer->GetProperty()->GetScalarSets();
    for ( int i = 0; i < (int)svs.size(); i++ )
    {
      ui->comboBoxScalarMap->addItem( svs[i].strName );
      if ( layer->GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarSet &&
           layer->GetProperty()->GetScalarSet() == i )
      {
        nSel = i + layers.size();
      }
    }
    ui->comboBoxScalarMap->addItem( "Load..." );
    if ( nSel >= 0 )
    {
      ui->comboBoxScalarMap->setCurrentIndex( nSel );
    }

    bShowSpline = layer->GetProperty()->GetShowSpline();
    ui->checkBoxShowSpline->setChecked( bShowSpline );
    ui->checkBoxSnapToCenter->setChecked( layer->GetProperty()->GetSnapToVoxelCenter() );
  }


  // MainWindow* mainWnd = MainWindow::GetMainWindowPointer();
  ui->colorpickerPointColor->setEnabled( layer );
  ui->comboBoxSplineColor->setEnabled( layer );

  ShowWidgets( m_widgetlistSpline, bShowSpline );
  ShowWidgets( m_widgetlistSolidColor, bShowSpline && layer && nColorMap == LayerPropertyPointSet::SolidColor );
  ShowWidgets( m_widgetlistHeatScale, bShowSpline && layer && nColorMap == LayerPropertyPointSet::HeatScale );

  BlockAllSignals( false );
}


void PanelPointSet::OnSliderOpacity( int nVal )
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    layer->GetProperty()->SetOpacity( nVal / 100.0 );
  }
}

void PanelPointSet::OnSliderMin(int nVal)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMin( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }

}

void PanelPointSet::OnSliderMid(int nVal)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMid( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelPointSet::OnSliderMax(int nVal)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    double fMin = layer->GetProperty()->GetScalarMinValue();
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleMax( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelPointSet::OnSliderOffset(int nVal)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    double fMax = layer->GetProperty()->GetScalarMaxValue();
    layer->GetProperty()->SetHeatScaleOffset( nVal / 100.0 * ( fMax + fMax ) - fMax );
  }
}

void PanelPointSet::OnLineEditMin(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && layer->GetProperty()->GetHeatScaleMin() != dVal )
  {
    layer->GetProperty()->SetHeatScaleMin( dVal );
  }
}

void PanelPointSet::OnLineEditMid(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && layer->GetProperty()->GetHeatScaleMid() != dVal )
  {
    layer->GetProperty()->SetHeatScaleMid( dVal );
  }
}

void PanelPointSet::OnLineEditMax(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && layer->GetProperty()->GetHeatScaleMax() != dVal )
  {
    layer->GetProperty()->SetHeatScaleMax( dVal );
  }
}

void PanelPointSet::OnLineEditOffset(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && layer->GetProperty()->GetHeatScaleOffset() != dVal )
  {
    layer->GetProperty()->SetHeatScaleOffset( dVal );
  }
}

void PanelPointSet::OnLineEditRadius(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && dVal > 0 && layer->GetProperty()->GetRadius() != dVal )
  {
    layer->GetProperty()->SetRadius( dVal );
  }
}

void PanelPointSet::OnLineEditSplineRadius(const QString& text)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  bool bOK;
  double dVal = text.toDouble( &bOK );
  if ( layer && bOK && dVal > 0 && layer->GetProperty()->GetSplineRadius() != dVal )
  {
    layer->GetProperty()->SetSplineRadius( dVal );
  }
}

void PanelPointSet::OnComboScalarMap(int nSel)
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(ui->comboBoxScalarMap->itemData( nSel ).value<QObject*>());
    if ( mri )
    {
      layer->GetProperty()->SetScalarLayer( mri );
    }
    else if ( nSel == ui->comboBoxScalarMap->count() - 1 )
    {
      LoadScalarValues();
    }
    else
    {
      int offset = ui->comboBoxScalarMap->count()-1 - layer->GetProperty()->GetNumberOfScalarSets();
      layer->GetProperty()->SetScalarSet( nSel - offset );
    }
  }
}

void PanelPointSet::LoadScalarValues()
{
  LayerPointSet* layer = GetCurrentLayer<LayerPointSet*>();
  if ( layer )
  {
    QString fn = QFileDialog::getOpenFileName( this, "Select scalar file",
                 "",
                 "All files (*)");
    if ( !fn.isEmpty() )
    {
      if ( !layer->GetProperty()->LoadScalarsFromFile( fn ) )
      {
        cout << "Load scalar values failed.\n";
      }
    }
    UpdateWidgets();
  }
}
