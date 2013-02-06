/**
 * @file  PanelVolume.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/02/06 18:35:43 $
 *    $Revision: 1.83 $
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
#include "PanelVolume.h"
#include "ui_PanelVolume.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerDTI.h"
#include "LayerPropertyDTI.h"
#include "LayerVolumeTrack.h"
#include "LUTDataHolder.h"
#include "MyUtils.h"
#include <QToolBar>
#include <QDebug>
#include <QStringList>
#include <QFileDialog>
#include <QMessageBox>
#include <QClipboard>
#include <QMimeData>

#define FS_VOLUME_SETTING_ID    "freesurfer/volume-setting"

PanelVolume::PanelVolume(QWidget *parent) :
  PanelLayer(parent),
  ui(new Ui::PanelVolume),
  m_curCTAB( NULL ),
  m_bShowExistingLabelsOnly(false)
{
  ui->setupUi(this);

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if ( !mainwnd )
  {
    return;
  }

  ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionNewVolume);
  ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionLoadVolume);
  ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionCloseVolume);
  ui->toolbar->insertAction(ui->actionMoveLayerUp, mainwnd->ui->actionSaveVolume);
  ui->toolbar->insertSeparator(ui->actionMoveLayerUp);

  m_luts = mainwnd->GetLUTData();

  m_widgetlistGrayScale << ui->checkBoxClearBackground
                        << ui->labelWindow
                        << ui->labelLevel
                        << ui->lineEditWindow
                        << ui->lineEditLevel
                        << ui->sliderWindow
                        << ui->sliderLevel;

  m_widgetlistHeatScale << ui->sliderMid
                        << ui->sliderOffset
                        << ui->lineEditMid
                        << ui->lineEditOffset
                        << ui->checkBoxClearHigher
                        << ui->checkBoxTruncate
                        << ui->checkBoxInvert
                        << ui->labelMid
                        << ui->labelOffset;

  m_widgetlistGenericColorMap << ui->lineEditMin
                              << ui->lineEditMax
                              << ui->sliderMin
                              << ui->sliderMax
                              << ui->labelMin
                              << ui->labelMax;

  m_widgetlistLUT << ui->treeWidgetColorTable
                  << ui->labelLookUpTable
                  << ui->comboBoxLookUpTable
                  << ui->colorLabelBrushValue
                  << ui->checkBoxShowExistingLabels;

  m_widgetlistDirectionCode << ui->comboBoxDirectionCode
                            << ui->labelDirectionCode;

  m_widgetlistFrame << ui->sliderFrame
                    << ui->spinBoxFrame
                    << ui->labelFrame
                    << ui->checkBoxRememberFrame
                    << ui->labelRememberFrame;

  m_widgetlistVector << ui->labelInversion
                     << ui->comboBoxInversion
                     << ui->labelRenderObject
                     << ui->comboBoxRenderObject
                     << ui->labelMask
                     << ui->comboBoxMask;

  m_widgetlistContour << ui->sliderContourThresholdLow
                      << ui->sliderContourThresholdHigh
                      << ui->lineEditContourThresholdLow
                      << ui->lineEditContourThresholdHigh
                      << ui->checkBoxUseColorMap
                      << ui->checkBoxContourExtractAll
                      << ui->colorPickerContour
                      << ui->labelContourThresholdHigh
                      << ui->labelContourThresholdLow
                      << ui->labelContourColor
                      << ui->sliderContourSmoothIteration
                      << ui->lineEditContourSmoothIteration
                      << ui->labelSmoothIteration
                      << ui->pushButtonContourSave
                      << ui->labelContourLabelRange
                      << ui->lineEditContourLabelRangeLow
                      << ui->lineEditContourLabelRangeHigh
                      << ui->checkBoxShowLabelContour
                      << ui->checkBoxUpsampleContour;

  m_widgetlistContourNormal << ui->sliderContourThresholdLow
      << ui->sliderContourThresholdHigh
      << ui->lineEditContourThresholdLow
      << ui->lineEditContourThresholdHigh
      << ui->checkBoxUseColorMap
      << ui->checkBoxContourExtractAll
      << ui->colorPickerContour
      << ui->labelContourThresholdHigh
      << ui->labelContourThresholdLow
      << ui->labelContourColor
      << ui->pushButtonContourSave;

  m_widgetlistContourLabel << ui->labelContourLabelRange
      << ui->lineEditContourLabelRangeLow
      << ui->lineEditContourLabelRangeHigh;

  m_widgetlistEditable << ui->labelBrushValue
                       << ui->lineEditBrushValue;

  m_widgetlistNormalDisplay << ui->labelOpacity
                            << ui->sliderOpacity
                            << ui->doubleSpinBoxOpacity
                            << ui->checkBoxSmooth
                            << ui->checkBoxUpsample
                            << ui->labelColorMap
                            << ui->comboBoxColorMap;

  m_widgetlistVolumeTrack << ui->treeWidgetColorTable << m_widgetlistFrame
                          << ui->labelSmoothIteration << ui->sliderContourSmoothIteration
                          << ui->lineEditContourSmoothIteration;

  m_widgetlistVolumeTrackSpecs
                        << ui->labelTrackVolumeThreshold
                        << ui->sliderTrackVolumeThresholdLow
                        << ui->lineEditTrackVolumeThresholdLow;

  QList<QWidget*> combo;
  combo << m_widgetlistGrayScale << m_widgetlistHeatScale
      << m_widgetlistGenericColorMap << m_widgetlistLUT
      << m_widgetlistDirectionCode << m_widgetlistVector
      << m_widgetlistContour << m_widgetlistEditable
      << ui->checkBoxSmooth << ui->checkBoxUpsample
      << ui->labelColorMap << ui->comboBoxColorMap
      << ui->checkBoxShowContour << ui->checkBoxShowOutline;

  combo = combo.toSet().toList();
  foreach (QWidget* w, m_widgetlistVolumeTrack)
  {
    int n = combo.indexOf(w);
    if (n >= 0)
      combo.removeAt(n);
  }
  m_widgetlistNonVolumeTrack = combo;

  LayerCollection* lc = mainwnd->GetLayerCollection("MRI");
  PanelLayer::InitializeLayerList( ui->treeWidgetLayers, lc );
  connect( ui->actionLockLayer, SIGNAL(toggled(bool)), lc, SLOT(LockCurrent(bool)) );
  connect( ui->actionMoveLayerUp, SIGNAL(triggered()), lc, SLOT(MoveLayerUp()));
  connect( ui->actionMoveLayerDown, SIGNAL(triggered()), lc, SLOT(MoveLayerDown()));
}

PanelVolume::~PanelVolume()
{
  delete ui;
}

void PanelVolume::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerMRI* layer = qobject_cast<LayerMRI*>(layer_in);
  if ( !layer )
  {
    return;
  }

  LayerPropertyMRI* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
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
  connect( ui->checkBoxClearBackground, SIGNAL(toggled(bool)), p, SLOT(SetClearZero(bool)) );
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
}

void PanelVolume::DoIdle()
{
  // update action status
  BlockAllSignals( true );
  QTreeWidgetItem* item = ui->treeWidgetLayers->currentItem();
  LayerMRI* layer = NULL;
  if ( item )
  {
    layer = qobject_cast<LayerMRI*>( item->data(0, Qt::UserRole).value<QObject*>() );
  }
  int nItemIndex = ui->treeWidgetLayers->indexOfTopLevelItem(item);
  ui->actionMoveLayerUp->setEnabled( item && !layer->IsLocked() && ui->treeWidgetLayers->topLevelItemCount() > 1 &&
                                     nItemIndex != 0 );
  ui->actionMoveLayerDown->setEnabled( item && !layer->IsLocked() && ui->treeWidgetLayers->topLevelItemCount() > 1 &&
                                       nItemIndex < ui->treeWidgetLayers->topLevelItemCount()-1 );
  ui->actionLockLayer->setEnabled( item );
  ui->actionLockLayer->setChecked( layer && layer->IsLocked() );
  ui->actionCopySetting->setEnabled(layer);
  QStringList strgs = qApp->clipboard()->text().split(",");
  bool bDataAvail = (!strgs.isEmpty() && strgs[0] == FS_VOLUME_SETTING_ID);
  ui->actionPasteSetting->setEnabled(layer && bDataAvail);
  ui->actionPasteSettingToAll->setEnabled(layer && bDataAvail);
  BlockAllSignals( false );
}

void PanelVolume::DoUpdateWidgets()
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

  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar &&
         allWidgets[i] != ui->toolbar2 && allWidgets[i]->parentWidget() != ui->toolbar2 )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  int nColorMap = LayerPropertyMRI::NoColorMap;
  ui->lineEditFileName->clear();
  if ( layer )
  {
    nColorMap = layer->GetProperty()->GetColorMap();
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );
    ui->checkBoxClearBackground->setChecked( layer->GetProperty()->GetClearZero() );
    if ( layer->IsTypeOf( "DTI" ) )
    {
      ui->lineEditFileName->setText(MyUtils::Win32PathProof(((LayerDTI*)layer)->GetVectorFileName()) );
    }
    else
    {
      ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
    }
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );

    ui->checkBoxSmooth->setChecked( layer->GetProperty()->GetTextureSmoothing() );
    ui->checkBoxUpsample->setChecked( layer->GetProperty()->GetUpSampleMethod() != 0 );
    ui->checkBoxUpsample->setEnabled( layer->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT );

    // color map settings
    ui->comboBoxColorMap->setCurrentIndex( layer->GetProperty()->GetColorMap() );
    ui->comboBoxLookUpTable->clear();
    for ( int i = 0; i < m_luts->GetCount(); i++ )
    {
      ui->comboBoxLookUpTable->addItem( m_luts->GetName( i ) );
    }
    if ( layer->GetEmbeddedColorTable() )
    {
      ui->comboBoxLookUpTable->addItem( "Embedded" );
    }
    ui->comboBoxLookUpTable->addItem( "Load lookup table..." );
    int nSel = m_luts->GetIndex( layer->GetProperty()->GetLUTCTAB() );
    ui->comboBoxLookUpTable->setCurrentIndex( nSel >= 0 ? nSel : m_luts->GetCount() );

    ChangeLineEditNumber( ui->lineEditBrushValue, layer->GetFillValue() );
    double dwindow = layer->GetProperty()->GetWindow();
    double dlevel  = layer->GetProperty()->GetLevel();
    ChangeLineEditNumber( ui->lineEditWindow, dwindow );
    ChangeLineEditNumber( ui->lineEditLevel, dlevel );
    double dMinTh = dlevel - dwindow/2;
    double dMaxTh = dlevel + dwindow/2;
    double dminvalue = layer->GetProperty()->GetMinValue();
    double dmaxvalue = layer->GetProperty()->GetMaxValue();
    double range_min = dminvalue - (dmaxvalue-dminvalue)/4;
    double range_max = dmaxvalue + (dmaxvalue-dminvalue)/4;
    if ( nColorMap == LayerPropertyMRI::Heat )
    {
      dMinTh = layer->GetProperty()->GetHeatScaleMinThreshold();
      dMaxTh = layer->GetProperty()->GetHeatScaleMaxThreshold();
      range_min = dminvalue;
      range_max = dmaxvalue;
    }
    else if ( nColorMap != LayerPropertyMRI::Grayscale && nColorMap != LayerPropertyMRI::LUT)
    {
      dMinTh = layer->GetProperty()->GetMinGenericThreshold();
      dMaxTh = layer->GetProperty()->GetMaxGenericThreshold();
    }
    double* windowrange = layer->GetProperty()->GetWindowRange();
    double* levelrange = layer->GetProperty()->GetLevelRange();
    ui->sliderWindow->setValue( (int)( ( dwindow - windowrange[0] ) / ( windowrange[1] - windowrange[0] ) * 100 ) );
    ui->sliderLevel->setValue( (int)( ( dlevel - levelrange[0] ) / ( levelrange[1] - levelrange[0] ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditMax, dMaxTh );
    ChangeLineEditNumber( ui->lineEditMin, dMinTh );
    ui->sliderMin->setValue( (int)( ( dMinTh - range_min ) / ( range_max - range_min ) * 100 ) );
    ui->sliderMax->setValue( (int)( ( dMaxTh - range_min ) / ( range_max - range_min ) * 100 ) );

    ui->sliderMid->setValue( (int)( ( layer->GetProperty()->GetHeatScaleMidThreshold() - range_min ) /
                                    ( range_max - range_min ) * 100 ) );
    ui->sliderOffset->setValue( (int)( ( layer->GetProperty()->GetHeatScaleOffset() + range_max ) /
                                       ( range_max + range_max ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditMid, layer->GetProperty()->GetHeatScaleMidThreshold() );
    ChangeLineEditNumber( ui->lineEditOffset, layer->GetProperty()->GetHeatScaleOffset() );
    ui->checkBoxClearHigher->setChecked( layer->GetProperty()->GetHeatScaleClearHigh() );
    ui->checkBoxTruncate->setChecked( layer->GetProperty()->GetHeatScaleTruncate() );
    ui->checkBoxInvert->setChecked( layer->GetProperty()->GetHeatScaleInvert() );

    ui->comboBoxColorMap->clear();
    ui->comboBoxColorMap->addItem( "Grayscale", LayerPropertyMRI::Grayscale );
    if ( layer->IsTypeOf( "DTI" ) )
    {
      ui->comboBoxColorMap->addItem( "Direction-coded", LayerPropertyMRI::DirectionCoded );
      ui->comboBoxDirectionCode->setCurrentIndex( ((LayerDTI*)layer)->GetProperty()->GetDirectionCode() );
    }
    else
    {
      ui->comboBoxColorMap->addItem( "Lookup Table", LayerPropertyMRI::LUT );
    }
    ui->comboBoxColorMap->addItem( "Heat", LayerPropertyMRI::Heat );
    ui->comboBoxColorMap->addItem( "Jet", LayerPropertyMRI::Jet );
    ui->comboBoxColorMap->addItem( "GE Color", LayerPropertyMRI::GEColor );
    ui->comboBoxColorMap->addItem( "NIH", LayerPropertyMRI::NIH );
    for ( int i = 0; i < ui->comboBoxColorMap->count(); i++ )
    {
      if ( ui->comboBoxColorMap->itemData( i ).toInt() == nColorMap )
      {
        ui->comboBoxColorMap->setCurrentIndex( i );
        break;
      }
    }

    int nFrames = layer->GetNumberOfFrames();
    if ( nFrames > 1 )
    {
      ui->sliderFrame->setRange( 1, nFrames );
      ui->spinBoxFrame->setRange( 1, nFrames );
      ui->checkBoxRememberFrame->setChecked(layer->GetProperty()->GetRememberFrameSettings());
    }
    ui->sliderFrame->setValue( layer->GetActiveFrame() + 1 );
    ChangeSpinBoxValue( ui->spinBoxFrame, layer->GetActiveFrame() + 1 );

    ui->checkBoxShowContour->setChecked( layer->GetProperty()->GetShowAsContour() );
    ui->sliderContourThresholdLow->setValue( (int)( ( layer->GetProperty()->GetContourMinThreshold() - dminvalue ) / ( dmaxvalue - dminvalue ) * 100 ) );
    ui->sliderContourThresholdHigh->setValue( (int)( ( layer->GetProperty()->GetContourMaxThreshold() - dminvalue ) / ( dmaxvalue - dminvalue ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditContourThresholdLow, layer->GetProperty()->GetContourMinThreshold() );
    ChangeLineEditNumber( ui->lineEditContourThresholdHigh, layer->GetProperty()->GetContourMaxThreshold() );
    ui->checkBoxUseColorMap->setChecked( layer->GetProperty()->GetContourUseImageColorMap() );
    ui->checkBoxUpsampleContour->setChecked( layer->GetProperty()->GetContourUpsample());
    ui->checkBoxContourExtractAll->setChecked( layer->GetProperty()->GetContourExtractAllRegions() );
    ui->sliderContourSmoothIteration->setValue( layer->GetProperty()->GetContourSmoothIterations() );
    ChangeLineEditNumber( ui->lineEditContourSmoothIteration, layer->GetProperty()->GetContourSmoothIterations() );
    double dTh1, dTh2;
    layer->GetProperty()->GetLabelContourRange(&dTh1, &dTh2);
    ChangeLineEditNumber( ui->lineEditContourLabelRangeLow, dTh1 );
    ChangeLineEditNumber( ui->lineEditContourLabelRangeHigh, dTh2 );

    ui->colorPickerContour->setEnabled( !layer->GetProperty()->GetContourUseImageColorMap() );
    double rgb[3];
    layer->GetProperty()->GetContourColor( rgb );
    ui->colorPickerContour->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );

    ui->comboBoxRenderObject->clear();
    if ( layer->GetProperty()->GetDisplayVector() )
    {
      ui->comboBoxRenderObject->addItem( "Simple Line" );
      ui->comboBoxRenderObject->addItem( "3D Bar (slow!)" );
      ui->comboBoxRenderObject->setCurrentIndex( layer->GetProperty()->GetVectorRepresentation() );
      ui->comboBoxInversion->setCurrentIndex( layer->GetProperty()->GetVectorInversion() );
    }
    else if ( layer->GetProperty()->GetDisplayTensor() )
    {
      ui->comboBoxRenderObject->addItem( "Boxoid" );
      ui->comboBoxRenderObject->addItem( "Ellipsoid (Very slow!)" );
      ui->comboBoxRenderObject->setCurrentIndex( layer->GetProperty()->GetTensorRepresentation() );
      ui->comboBoxInversion->setCurrentIndex( layer->GetProperty()->GetTensorInversion() );
    }

    ui->checkBoxShowInfo->setChecked( layer->GetProperty()->GetShowInfo() );
    ui->checkBoxProjectionMap->setChecked( layer->GetProperty()->GetShowProjectionMap());

    ui->checkBoxShowOutline->setChecked( layer->GetProperty()->GetShowLabelOutline() );
    ui->checkBoxShowOutline->setVisible( !layer->IsTypeOf("DTI") );
    ui->checkBoxProjectionMap->setVisible(!layer->IsTypeOf("DTI") );

    //    ui->m_choiceUpSampleMethod->SetSelection( layer->GetProperty()->GetUpSampleMethod() );
    ui->checkBoxShowExistingLabels->setEnabled(!layer->GetAvailableLabels().isEmpty());
  }

  bool bNormalDisplay = (layer && !layer->GetProperty()->GetDisplayVector() && !layer->GetProperty()->GetDisplayTensor());

  if (layer && layer->IsTypeOf("VolumeTrack"))
  {
    ShowWidgets(m_widgetlistNonVolumeTrack, false);
    ShowWidgets(m_widgetlistVolumeTrack, true);
    if (m_curCTAB != layer->GetEmbeddedColorTable())
      PopulateColorTable( layer->GetEmbeddedColorTable() );
  }
  else
  {
    ShowWidgets( m_widgetlistNormalDisplay, bNormalDisplay );
    ShowWidgets( m_widgetlistGrayScale, bNormalDisplay && nColorMap == LayerPropertyMRI::Grayscale );
    ShowWidgets( m_widgetlistHeatScale, bNormalDisplay && nColorMap == LayerPropertyMRI::Heat );
    ShowWidgets( m_widgetlistGenericColorMap, bNormalDisplay && nColorMap != LayerPropertyMRI::LUT &&
                 nColorMap != LayerPropertyMRI::DirectionCoded );
    ShowWidgets( m_widgetlistLUT, bNormalDisplay && nColorMap == LayerPropertyMRI::LUT );
    ShowWidgets( m_widgetlistDirectionCode, bNormalDisplay && nColorMap == LayerPropertyMRI::DirectionCoded );
    ShowWidgets( m_widgetlistEditable, bNormalDisplay && layer->IsEditable() );
    ShowWidgets( m_widgetlistFrame, layer &&
                 !layer->IsTypeOf( "DTI" ) &&
                 layer->GetNumberOfFrames() > 1 );

    ui->sliderFrame->setEnabled( layer &&
                                 !layer->GetProperty()->GetDisplayVector() &&
                                 !layer->GetProperty()->GetDisplayTensor() );
    ui->spinBoxFrame->setEnabled( layer &&
                                  !layer->GetProperty()->GetDisplayVector() &&
                                  !layer->GetProperty()->GetDisplayTensor() );
    ui->checkBoxDisplayVector->setVisible( layer && ( layer->IsTypeOf( "DTI" ) || layer->GetNumberOfFrames() == 3 ) );
    ui->checkBoxDisplayVector->setChecked( layer && layer->GetProperty()->GetDisplayVector() );
    ui->checkBoxDisplayTensor->setVisible( layer && layer->GetNumberOfFrames() == 9 );
    ui->checkBoxDisplayTensor->setChecked( layer && layer->GetProperty()->GetDisplayTensor() );
    ShowWidgets( m_widgetlistVector, ui->checkBoxDisplayVector->isChecked() || ui->checkBoxDisplayTensor->isChecked() );
    ShowWidgets( m_widgetlistContour, ui->checkBoxShowContour->isChecked() );
    ui->checkBoxShowContour->setVisible( bNormalDisplay );
    if (layer && ui->checkBoxShowContour->isChecked())
    {
      ui->checkBoxShowLabelContour->setChecked(layer->GetProperty()->GetShowAsLabelContour());
      ShowWidgets( m_widgetlistContourNormal, !layer->GetProperty()->GetShowAsLabelContour());
      ShowWidgets( m_widgetlistContourLabel, layer->GetProperty()->GetShowAsLabelContour());
    }

    //  ShowWidgets( m_widgetlistContour, false );
    //  m_checkContour->Show( false /*nColorMap == LayerPropertyMRI::LUT*/ );

    if ( layer && layer->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT )
    {
      if ( m_curCTAB != layer->GetProperty()->GetLUTCTAB() ||
           m_bShowExistingLabelsOnly != ui->checkBoxShowExistingLabels->isChecked())
      {
        PopulateColorTable( layer->GetProperty()->GetLUTCTAB() );
      }

      for ( int i = 0; i < ui->treeWidgetColorTable->topLevelItemCount(); i++ )
      {
        QTreeWidgetItem* item = ui->treeWidgetColorTable->topLevelItem( i );
        QStringList strglist = item->text(0).split( " " );
        bool bOK;
        double dvalue = strglist[0].trimmed().toDouble( &bOK );
        if ( bOK && dvalue == layer->GetFillValue() )
        {
          ui->treeWidgetColorTable->setCurrentItem( item );
          break;
        }
      }
      UpdateColorLabel();
    }
  }

  UpdateTrackVolumeThreshold();

  ui->checkBoxUpsampleContour->hide();

  BlockAllSignals( false );
}

void PanelVolume::OnColorTableCurrentItemChanged( QTreeWidgetItem* item )
{
  if ( item )
  {
    QStringList strglist = item->text( 0 ).split(" ");
    double val = strglist[0].toDouble();
    LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
    if ( layer->IsTypeOf("VolumeTrack") )
    {
      UpdateTrackVolumeThreshold();
    }
    else
    {
      layer->SetFillValue( val );
    }
    ChangeLineEditNumber( ui->lineEditBrushValue, val );
    UpdateColorLabel();
  }
}

void PanelVolume::UpdateColorLabel()
{
  QTreeWidgetItem* item = ui->treeWidgetColorTable->currentItem();
  if ( item )
  {
    QColor color = item->data( 0, Qt::UserRole ).value<QColor>();
    if ( color.isValid() )
    {
      QPixmap pix( 30, 20 );
      pix.fill( color );
      ui->colorLabelBrushValue->setPixmap( pix );
    }
  }
  else
  {
    ui->colorLabelBrushValue->setPixmap( QPixmap() );
  }
}

void PanelVolume::UpdateTrackVolumeThreshold()
{
  LayerVolumeTrack* layer = GetCurrentLayer<LayerVolumeTrack*>();
  QTreeWidgetItem* item = ui->treeWidgetColorTable->currentItem();
  if ( item && layer )
  {
    int nLabel = item->text(0).split(" ").at(0).toInt();
    ui->sliderTrackVolumeThresholdLow->blockSignals(true);
    ui->lineEditTrackVolumeThresholdLow->blockSignals(true);
    double fMin = layer->GetProperty()->GetMinValue();
    double fMax = layer->GetProperty()->GetMaxValue()/4;
    ui->sliderTrackVolumeThresholdLow->setValue( (int)( ( layer->GetThreshold(nLabel) - fMin ) / ( fMax - fMin ) * 100 ) );
    ChangeLineEditNumber( ui->lineEditTrackVolumeThresholdLow, layer->GetThreshold(nLabel) );
    ui->sliderTrackVolumeThresholdLow->blockSignals(false);
    ui->lineEditTrackVolumeThresholdLow->blockSignals(false);
    layer->Highlight(nLabel);
  }
  ShowWidgets(m_widgetlistVolumeTrackSpecs, layer);
  EnableWidgets(this->m_widgetlistVolumeTrackSpecs, item);
}

void PanelVolume::PopulateColorTable( COLOR_TABLE* ct )
{
  if ( ct && ct != m_curCTAB )
  {
    m_curCTAB = ct;
    ui->treeWidgetColorTable->clear();
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    int nValid = 0;
    char name[1000];
    int nSel = -1;
    int nValue = 0;
    bool bOK;
    nValue = ui->lineEditBrushValue->text().toInt( &bOK );
    LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
    if ( !bOK && layer )
    {
      nValue = (int)layer->GetFillValue();
    }

    QList<int> labels;
    if (layer)
      labels = layer->GetAvailableLabels();
    int nValidCount = 0;
    for ( int i = 0; i < nTotalCount; i++ )
    {
      CTABisEntryValid( ct, i, &nValid );
      if ( nValid )
      {
        CTABcopyName( ct, i, name, 1000 );
        QTreeWidgetItem* item = new QTreeWidgetItem( ui->treeWidgetColorTable );
        item->setText( 0, QString("%1 %2").arg(i).arg(name) );
        item->setToolTip( 0, name );
        int nr, ng, nb;
        CTABrgbAtIndexi( ct, i, &nr, &ng, &nb );
        QColor color( nr, ng, nb );
        QPixmap pix(13, 13);
        pix.fill( color );
        item->setIcon(0, QIcon(pix) );
        item->setData( 0, Qt::UserRole, color );
        item->setData(0, Qt::UserRole+1, i);
        if (m_bShowExistingLabelsOnly && !labels.isEmpty())
        {
          item->setHidden(!labels.contains(i));
        }
        if ( i == nValue )
        {
          nSel = nValidCount;
        }
        nValidCount++;
      }
    }
    if ( nSel >= 0 )
    {
      ui->treeWidgetColorTable->setCurrentItem( ui->treeWidgetColorTable->topLevelItem( nSel ) );
    }
  }
}

void PanelVolume::OnLineEditBrushValue( const QString& strg )
{
  QString text = strg.trimmed();
  bool bOK;
  int nVal = text.toInt( &bOK );
  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  QList<int> labels = layer->GetAvailableLabels();
  if ( text.isEmpty() )
  {
    for ( int i = 0; i < ui->treeWidgetColorTable->topLevelItemCount(); i++ )
    {
      QTreeWidgetItem* item = ui->treeWidgetColorTable->topLevelItem( i );
      if (m_bShowExistingLabelsOnly)
      {
        int n = item->data(0, Qt::UserRole+1).toInt();
        item->setHidden(!labels.contains(n));
      }
      else
        item->setHidden( false );
    }
  }
  else if ( bOK )
  {
    if ( layer )
    {
      layer->SetFillValue( nVal );
    }
    bool bFound = false;
    for ( int i = 0; i < ui->treeWidgetColorTable->topLevelItemCount(); i++ )
    {
      QTreeWidgetItem* item = ui->treeWidgetColorTable->topLevelItem( i );
      if (m_bShowExistingLabelsOnly)
      {
        int n = item->data(0, Qt::UserRole+1).toInt();
        item->setHidden(!labels.contains(n));
      }
      else
        item->setHidden( false );
      QStringList strglist = item->text(0).split( " " );
      if ( strglist[0].toDouble() == layer->GetFillValue() )
      {
        ui->treeWidgetColorTable->setCurrentItem( item );
        bFound = true;
      }
    }
    if( !bFound )
    {
      ui->treeWidgetColorTable->setCurrentItem( NULL );
    }
    UpdateColorLabel();
  }
  else
  {
    for ( int i = 0; i < ui->treeWidgetColorTable->topLevelItemCount(); i++ )
    {
      QTreeWidgetItem* item = ui->treeWidgetColorTable->topLevelItem( i );
      if ( item->text(0).contains( text, Qt::CaseInsensitive ) )
      {
        if (m_bShowExistingLabelsOnly)
        {
          int n = item->data(0, Qt::UserRole+1).toInt();
          item->setHidden(!labels.contains(n));
        }
        else
          item->setHidden( false );
      }
      else
      {
        item->setHidden( true );
      }
    }
  }
}

void PanelVolume::OnComboColorMap( int nSel )
{
  /*
  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  if ( layer && nSel >= 0 )
  {
    nSel = ui->comboBoxColorMap->itemData(nSel).toInt();
    layer->GetProperty()->SetColorMap( nSel );
  }
  */
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    if (nSel >= 0)
    {
      nSel = ui->comboBoxColorMap->itemData(nSel).toInt();
      layer->GetProperty()->SetColorMap( nSel );
    }
  }
}

void PanelVolume::OnComboLookupTable( int nSel )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    if ( nSel == ui->comboBoxLookUpTable->count()-1 )
    {
      MainWindow::GetMainWindow()->LoadLUT();
    }
    else
    {
      if ( nSel < m_luts->GetCount() )
      {
        COLOR_TABLE* ct = m_luts->GetColorTable( nSel );
        layer->GetProperty()->SetLUTCTAB( ct );
      }
      else
      {
        layer->GetProperty()->SetLUTCTAB( layer->GetEmbeddedColorTable() );
      }
    }
  }
}

void PanelVolume::OnCheckShowContour(bool bShow)
{
  ShowWidgets( m_widgetlistContour, bShow );
}

void PanelVolume::OnCheckShowLabelContour(bool bShow)
{
  ShowWidgets( m_widgetlistContourNormal, !bShow);
  ShowWidgets( m_widgetlistContourLabel, bShow);
}

void PanelVolume::OnSliderOpacity( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    layer->GetProperty()->SetOpacity( nVal/100.0 );
  }
}

void PanelVolume::OnSliderWindow( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double* r = layer->GetProperty()->GetWindowRange();
    layer->GetProperty()->SetWindow( nVal / 100.0 * ( r[1] - r[0] ) + r[0] );
  }
}

void PanelVolume::OnSliderLevel( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double* r = layer->GetProperty()->GetLevelRange();
    layer->GetProperty()->SetLevel( nVal / 100.0 * ( r[1] - r[0] ) + r[0] );
  }
}

void PanelVolume::OnSliderMin( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double fMin = layer->GetProperty()->GetMinValue();
    double fMax = layer->GetProperty()->GetMaxValue();
    double fScaleMin = fMin - (fMax-fMin)/4;
    double fScaleMax = fMax + (fMax-fMin)/4;
    switch ( layer->GetProperty()->GetColorMap() )
    {
    case LayerPropertyMRI::Grayscale:
      layer->GetProperty()->SetMinGrayscaleWindow( nVal /
          100.0 * ( fScaleMax - fScaleMin ) + fScaleMin );
      break;
    case LayerPropertyMRI::Heat:
      layer->GetProperty()->SetHeatScaleMinThreshold( nVal /
          100.0 * ( fMax - fMin ) + fMin );
      break;
    default:
      layer->GetProperty()->SetMinGenericThreshold( nVal /
          100.0 * ( fScaleMax - fScaleMin ) + fScaleMin  );
      break;
    }
  }
}

void PanelVolume::OnSliderMid( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double fMin = layer->GetProperty()->GetMinValue();
    double fMax = layer->GetProperty()->GetMaxValue();
    layer->GetProperty()->SetHeatScaleMidThreshold( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelVolume::OnSliderMax( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double fMin = layer->GetProperty()->GetMinValue();
    double fMax = layer->GetProperty()->GetMaxValue();
    double fScaleMin = fMin - (fMax-fMin)/4;
    double fScaleMax = fMax + (fMax-fMin)/4;
    switch ( layer->GetProperty()->GetColorMap() )
    {
    case LayerPropertyMRI::Grayscale:
      layer->GetProperty()->SetMaxGrayscaleWindow( nVal /
          100.0 * ( fScaleMax - fScaleMin ) + fScaleMin );
      break;
    case LayerPropertyMRI::Heat:
      layer->GetProperty()->SetHeatScaleMaxThreshold( nVal /
          100.0 * ( fMax - fMin ) + fMin );
      break;
    default:
      layer->GetProperty()->SetMaxGenericThreshold( nVal /
          100.0 * ( fScaleMax - fScaleMin ) + fScaleMin  );
      break;
    }
  }
}

void PanelVolume::OnSliderOffset( int nVal )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    double fMax = layer->GetProperty()->GetMaxValue();
    layer->GetProperty()->SetHeatScaleOffset( nVal / 100.0 * ( fMax + fMax ) - fMax );
  }
}

void PanelVolume::OnLineEditWindow( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if (bOK && layer->GetProperty()->GetWindow() != dVal )
    {
      layer->GetProperty()->SetWindow( dVal );
    }
  }
}

void PanelVolume::OnLineEditLevel( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( bOK && layer->GetProperty()->GetLevel() != dVal )
    {
      layer->GetProperty()->SetLevel( dVal );
    }
  }
}

void PanelVolume::OnLineEditMin( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK )
    {
      switch ( layer->GetProperty()->GetColorMap() )
      {
      case LayerPropertyMRI::Grayscale:
        layer->GetProperty()->SetMinGrayscaleWindow( dVal );
        break;
      case LayerPropertyMRI::Heat:
        layer->GetProperty()->SetHeatScaleMinThreshold( dVal );
        break;
      default:
        layer->GetProperty()->SetMinGenericThreshold( dVal );
        break;
      }
    }
  }
}

void PanelVolume::OnLineEditMid( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleMidThreshold() != dVal )
    {
      layer->GetProperty()->SetHeatScaleMidThreshold( dVal );
    }
  }
}

void PanelVolume::OnLineEditMax( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK )
    {
      switch ( layer->GetProperty()->GetColorMap() )
      {
      case LayerPropertyMRI::Grayscale:
        layer->GetProperty()->SetMaxGrayscaleWindow( dVal );
        break;
      case LayerPropertyMRI::Heat:
        layer->GetProperty()->SetHeatScaleMaxThreshold( dVal );
        break;
      default:
        layer->GetProperty()->SetMaxGenericThreshold( dVal );
        break;
      }
    }
  }
}

void PanelVolume::OnLineEditOffset( const QString& text )
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    bool bOK;
    double dVal = text.toDouble( &bOK );
    if ( layer && bOK && layer->GetProperty()->GetHeatScaleOffset() != dVal )
    {
      layer->GetProperty()->SetHeatScaleOffset( dVal );
    }
  }
}

void PanelVolume::OnSliderContourMin(int nval)
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    if ( layer )
    {
      double fMin = layer->GetProperty()->GetMinValue();
      double fMax = layer->GetProperty()->GetMaxValue();
      ChangeLineEditNumber( ui->lineEditContourThresholdLow,
                            nval / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourMax(int nval)
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    if ( layer )
    {
      double fMin = layer->GetProperty()->GetMinValue();
      double fMax = layer->GetProperty()->GetMaxValue();
      ChangeLineEditNumber( ui->lineEditContourThresholdHigh,
                            nval / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnSliderContourSmooth(int nval)
{
  ChangeLineEditNumber(ui->lineEditContourSmoothIteration, nval);
}

void PanelVolume::OnContourValueChanged()
{
  bool bOK;
  double fMin, fMax = 0;
  int nSmooth = 30;
  if (ui->checkBoxShowLabelContour->isChecked())
  {
    fMin = ui->lineEditContourLabelRangeLow->text().trimmed().toDouble(&bOK);
    if (bOK)
    {
      fMax = ui->lineEditContourLabelRangeHigh->text().trimmed().toDouble(&bOK);
    }
    if (bOK)
    {
      nSmooth = ui->lineEditContourSmoothIteration->text().trimmed().toInt(&bOK);
    }
    QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
    foreach (LayerMRI* layer, layers)
    {
      if (layer && bOK)
      {
        if (sender() == ui->lineEditContourSmoothIteration ||
            sender() == ui->sliderContourSmoothIteration )
        {
          layer->GetProperty()->SetContourSmoothIterations(nSmooth);
        }
        else
        {
          layer->GetProperty()->SetLabelContourRange(fMin, fMax);
        }
      }
    }
  }
  else
  {
    fMin = ui->lineEditContourThresholdLow->text().trimmed().toDouble(&bOK);
    if (bOK)
    {
      fMax = ui->lineEditContourThresholdHigh->text().trimmed().toDouble(&bOK);
    }
    if (bOK)
    {
      nSmooth = ui->lineEditContourSmoothIteration->text().trimmed().toInt(&bOK);
    }
    QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
    foreach (LayerMRI* layer, layers)
    {
      if (layer && bOK)
      {
        if (sender() == ui->lineEditContourSmoothIteration ||
            sender() == ui->sliderContourSmoothIteration )
        {
          layer->GetProperty()->SetContourSmoothIterations(nSmooth);
        }
        else
        {
          layer->GetProperty()->SetContourThreshold(fMin, fMax);
        }
      }
    }
  }
}

void PanelVolume::OnContourSave()
{
  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  if ( layer )
  {
    QString fn = QFileDialog::getSaveFileName( this,
                 "Save iso-surface",
                 MainWindow::GetMainWindow()->AutoSelectLastDir("mri") + "/" + layer->GetName() + ".vtk",
                 "VTK files (*.vtk);;All files (*)");
    if ( !fn.isEmpty() )
    {
      if ( !layer->SaveContourToFile( fn ) )
      {
        QMessageBox::warning(this, "Error", "Can not save surface to file.");
      }
    }
  }
}

void PanelVolume::OnSliderTrackVolumeMin(int nval)
{
  QList<LayerMRI*> layers = GetSelectedLayers<LayerMRI*>();
  foreach (LayerMRI* layer, layers)
  {
    if ( layer && layer->IsTypeOf("VolumeTrack"))
    {
      double fMin = layer->GetProperty()->GetMinValue();
      double fMax = layer->GetProperty()->GetMaxValue()/4;
      ChangeLineEditNumber( ui->lineEditTrackVolumeThresholdLow,
                            nval / 100.0 * ( fMax - fMin ) + fMin );
    }
  }
}

void PanelVolume::OnTrackVolumeThresholdChanged()
{
  QTreeWidgetItem* item = ui->treeWidgetColorTable->currentItem();
  if (!item)
    return;

  int nLabel = item->text(0).split(" ").at(0).toInt();
  bool bOK;
  double fMin;
  fMin = ui->lineEditTrackVolumeThresholdLow->text().trimmed().toDouble(&bOK);
  if (bOK)
  {
    LayerVolumeTrack* layer = GetCurrentLayer<LayerVolumeTrack*>();
    if (layer)
    {
      layer->SetThreshold(nLabel, fMin);
    }
  }
}

void PanelVolume::OnCopySettings()
{
  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  if ( layer )
  {
    QVariantMap map = layer->GetProperty()->GetActiveSettings();
    QStringList keys = map.keys();
    QStringList strgs;
    strgs << FS_VOLUME_SETTING_ID;
    for (int i = 0; i < keys.size(); i++ )
    {
      strgs << keys[i] << map[keys[i]].toString();
    }
    qApp->clipboard()->setText(strgs.join(","));
  }
}

void PanelVolume::OnPasteSettings()
{
  LayerMRI* layer = GetCurrentLayer<LayerMRI*>();
  if ( layer )
  {
    QStringList strgs = qApp->clipboard()->text().split(",");
    if (!strgs.isEmpty() && strgs[0] == FS_VOLUME_SETTING_ID)
    {
      QVariantMap map;
      for (int i = 1; i < strgs.size(); i+=2 )
      {
        bool bOK;
        double val = strgs[i+1].toDouble(&bOK);
        if (bOK)
        {
          map[strgs[i]] = val;
        }
      }
      layer->GetProperty()->RestoreSettings(map);
    }
  }
}

void PanelVolume::OnPasteSettingsToAll()
{
  QStringList strgs = qApp->clipboard()->text().split(",");
  if (!strgs.isEmpty() && strgs[0] == FS_VOLUME_SETTING_ID)
  {
    QVariantMap map;
    for (int i = 1; i < strgs.size(); i+=2 )
    {
      bool bOK;
      double val = strgs[i+1].toDouble(&bOK);
      if (bOK)
      {
        map[strgs[i]] = val;
      }
    }
    QList<Layer*> layers = GetSelectedLayers<Layer*>();
    if (layers.size() < 2)
      layers = MainWindow::GetMainWindow()->GetLayerCollection("MRI")->GetLayers();
    for (int i = 0; i < layers.size(); i++)
    {
      ((LayerMRI*)layers[i])->GetProperty()->RestoreSettings(map);
    }
  }
}

void PanelVolume::OnActiveFrameChanged(int nFrame)
{
  LayerVolumeTrack* layer = GetCurrentLayer<LayerVolumeTrack*>();
  if ( layer )
  {
    int nLabel = layer->GetFrameLabel(nFrame);
    for (int i = 0; i < ui->treeWidgetColorTable->topLevelItemCount(); i++)
    {
      QTreeWidgetItem* item = ui->treeWidgetColorTable->topLevelItem(i);
      if ( item->text(0).split(" ").at(0).toInt() == nLabel )
      {
        ui->treeWidgetColorTable->setCurrentItem(item);
        return;
      }
    }
  }
}

void PanelVolume::OnShowExistingLabelsOnly(bool b)
{
  m_bShowExistingLabelsOnly = b;
//  this->UpdateWidgets();
  OnLineEditBrushValue(ui->lineEditBrushValue->text());
}
