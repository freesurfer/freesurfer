/**
 * @file  PanelSurface.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2017/02/01 15:28:54 $
 *    $Revision: 1.74 $
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
#include "PanelSurface.h"
#include "MainWindow.h"
#include "ui_PanelSurface.h"
#include "ui_MainWindow.h"
#include <QToolBar>
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include "LayerCollection.h"
#include "FSSurface.h"
#include "SurfaceAnnotation.h"
#include "SurfaceOverlay.h"
#include "SurfaceLabel.h"
#include "SurfaceSpline.h"
#include "WindowConfigureOverlay.h"
#include "MyUtils.h"
#include <QMessageBox>
#include <QDebug>
#include "SurfacePath.h"
#include <QToolButton>
#include "LayerMRI.h"
#include "DialogCustomFill.h"
#include <QFileDialog>
#include "DialogSurfaceLabelOperations.h"

PanelSurface::PanelSurface(QWidget *parent) :
  PanelLayer("Surface", parent),
  ui(new Ui::PanelSurface)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  ui->toolbar->insertAction(ui->actionShowOverlay, mainwnd->ui->actionLoadSurface);
  ui->toolbar->insertAction(ui->actionShowOverlay, mainwnd->ui->actionCloseSurface);
  ui->toolbar->insertSeparator(ui->actionShowOverlay);
  m_toolButtonSurface = new QToolButton(this);
  m_toolButtonSurface->setStyleSheet("QToolButton::menu-indicator{image: none;}");
  m_toolButtonSurface->setPopupMode(QToolButton::InstantPopup);
  m_toolButtonSurface->setIcon(QIcon(":/resource/icons/surface_main.png"));
  ui->toolbar->insertWidget(ui->actionShowOverlay, m_toolButtonSurface);
  ui->toolbarSurfaces->hide();

  QMenu* menu = new QMenu(this);
  menu->addAction(ui->actionSurfaceMain);
  menu->addAction(ui->actionSurfaceInflated);
  menu->addAction(ui->actionSurfaceWhite);
  menu->addAction(ui->actionSurfacePial);
  menu->addAction(ui->actionSurfaceOriginal);
  m_toolButtonSurface->setMenu(menu);

  //  ui->treeWidgetLabels->hide();

  m_widgetsSlope << ui->sliderSlope
                 << ui->lineEditSlope
                 << ui->labelSlope;

  m_widgetsMidPoint << ui->sliderMidPoint
                    << ui->lineEditMidPoint
                    << ui->labelMidPoint;

  m_widgetsVector << ui->colorpickerVectorColor
                  << ui->spinBoxVectorPointSize
                  << ui->labelVectorColor
                  << ui->labelVectorPointSize;

  m_widgetsVertex << ui->colorpickerVertexColor
                  << ui->spinBoxVertexPointSize
                  << ui->labelVertexColor
                  << ui->labelVertexPointSize;

  m_widgetsMesh << ui->comboBoxMeshColor
                << ui->labelMeshColor;

  m_widgetsLabel << ui->colorpickerLabelColor
                 << ui->labelLabelColor
                 << ui->checkBoxLabelOutline
                 << ui->treeWidgetLabels
                 << ui->widgetLabelOrderButtons
                 << ui->labelLabelThreshold
                 << ui->comboBoxLabelColorCode
                 << ui->lineEditLabelThreshold
                 << ui->lineEditLabelHeatscaleMin
                 << ui->lineEditLabelHeatscaleMax
                 << ui->labelHeatscaleRange
                 << ui->labelLabelZOrderLabel
                 << ui->spinBoxZOrderLabel
                 << ui->labelLabelOpacity
                 << ui->lineEditLabelOpacity
                 << ui->toolButtonMoreLabelOptions;

  m_widgetsOverlay << ui->pushButtonConfigureOverlay
                   << ui->pushButtonRemoveOverlay
                   << ui->labelLabelZOrderOverlay
                   << ui->spinBoxZOrderOverlay;

  m_widgetsAnnotation << ui->checkBoxAnnotationOutline
                      << ui->labelLabelZOrderAnnotation
                      << ui->spinBoxZOrderAnnotation;

  m_widgetsSpline << ui->colorpickerSplineColor
                  << ui->labelSplineColor
                  << ui->checkBoxSplineProjection
                  << ui->treeWidgetSplines;

  ui->actionSurfaceMain->setData( FSSurface::SurfaceMain );
  ui->actionSurfaceInflated->setData( FSSurface::SurfaceInflated );
  ui->actionSurfaceOriginal->setData( FSSurface::SurfaceOriginal );
  ui->actionSurfacePial->setData( FSSurface::SurfacePial );
  ui->actionSurfaceWhite->setData( FSSurface::SurfaceWhite );
  QActionGroup* ag = new QActionGroup( this );
  ag->addAction(ui->actionSurfaceMain);
  ag->addAction(ui->actionSurfaceInflated);
  ag->addAction(ui->actionSurfaceWhite);
  ag->addAction(ui->actionSurfacePial);
  ag->addAction(ui->actionSurfaceOriginal);
  ag->setExclusive( true );
  connect( ag, SIGNAL(triggered(QAction*)), this, SLOT(OnChangeSurfaceType(QAction*)));
  m_actGroupSurface = ag;

  LayerCollection* lc = mainwnd->GetLayerCollection("Surface");
  connect( ui->actionLockLayer, SIGNAL(toggled(bool)), this, SLOT(OnLockLayer(bool)) );
  connect( ui->actionMoveLayerUp, SIGNAL(triggered()), lc, SLOT(MoveLayerUp()));
  connect( ui->actionMoveLayerDown, SIGNAL(triggered()), lc, SLOT(MoveLayerDown()));

  m_wndConfigureOverlay = new WindowConfigureOverlay( this );
  m_wndConfigureOverlay->hide();
  connect( mainwnd->GetLayerCollection("Surface"), SIGNAL(ActiveLayerChanged(Layer*)),
           m_wndConfigureOverlay, SLOT(OnActiveSurfaceChanged(Layer*)));
  connect(m_wndConfigureOverlay, SIGNAL(ActiveFrameChanged()), mainwnd, SLOT(UpdateInfoPanel()));
  connect(mainwnd, SIGNAL(SlicePositionChanged()), m_wndConfigureOverlay, SLOT(OnCurrentVertexChanged()));
  connect(m_wndConfigureOverlay, SIGNAL(MaskLoadRequested(QString)), mainwnd, SLOT(OnLoadSurfaceLabelRequested(QString)));
  connect(m_wndConfigureOverlay, SIGNAL(OverlayChanged()), SLOT(UpdateWidgets()));

  connect(ui->checkBoxLabelOutline, SIGNAL(toggled(bool)), this, SLOT(OnCheckBoxLabelOutline(bool)));
  connect(ui->colorpickerLabelColor, SIGNAL(colorChanged(QColor)), this, SLOT(OnColorPickerLabelColor(QColor)));
  connect(ui->treeWidgetLabels, SIGNAL(MenuGoToCentroid()), mainwnd, SLOT(OnGoToSurfaceLabel()));
  connect(ui->treeWidgetLabels, SIGNAL(MenuResample()), this, SLOT(OnLabelResample()));
  connect(ui->treeWidgetLabels, SIGNAL(MenuMoreOps()), this, SLOT(OnLabelMoreOps()));
  connect(ui->treeWidgetLabels, SIGNAL(MenuSaveAs()), this, SLOT(OnSaveLabelAs()));

  connect(ui->actionCut, SIGNAL(toggled(bool)), SLOT(OnButtonEditCut(bool)));
  connect(ui->actionPath, SIGNAL(toggled(bool)), SLOT(OnButtonEditPath(bool)));

  connect(ui->actionCutLine, SIGNAL(triggered(bool)), SLOT(OnButtonCutLine()));
  connect(ui->actionCutClosedLine, SIGNAL(triggered(bool)), SLOT(OnButtonCutClosedLine()));
  connect(ui->actionCutClear, SIGNAL(triggered(bool)), SLOT(OnButtonClearCuts()));
  connect(ui->actionFillUncutArea, SIGNAL(triggered(bool)), SLOT(OnButtonFillUncutArea()));
  connect(ui->actionUndoCut, SIGNAL(triggered(bool)), SLOT(OnButtonUndoCut()));

  connect(ui->actionMakePath, SIGNAL(triggered(bool)), SLOT(OnButtonMakePath()));
  connect(ui->actionMakeClosedPath, SIGNAL(triggered(bool)), SLOT(OnButtonMakeClosedPath()));
  connect(ui->actionDeletePath, SIGNAL(triggered(bool)), SLOT(OnButtonDeletePath()));
  connect(ui->actionPathFill, SIGNAL(triggered(bool)), SLOT(OnButtonCustomFillPath()));

  connect(ui->actionClearMarks, SIGNAL(triggered(bool)), SLOT(OnButtonClearMarks()));

  connect(ui->toolButtonLabelUp, SIGNAL(clicked(bool)), SLOT(OnButtonLabelUp()));
  connect(ui->toolButtonLabelDown, SIGNAL(clicked(bool)), SLOT(OnButtonLabelDown()));

  m_dlgCustomFill = new DialogCustomFill(this);
  m_dlgCustomFill->hide();
  connect(m_dlgCustomFill, SIGNAL(CustomFillTriggered(QVariantMap)), SLOT(OnCustomFillTriggered(QVariantMap)));

  m_dlgLabelOps = new DialogSurfaceLabelOperations(this);
  m_dlgLabelOps->hide();
  connect(m_dlgLabelOps, SIGNAL(OperationTriggered(QVariantMap)), SLOT(OnLabelOperation(QVariantMap)));

  ag = new QActionGroup(this);
  ag->addAction(ui->actionCutLine);
  ag->addAction(ui->actionCutClosedLine);
  ag->addAction(ui->actionFillUncutArea);
  ag->addAction(ui->actionCutClear);
  ag->addAction(ui->actionUndoCut);
  ag->setVisible(false);
  connect(ui->actionCut, SIGNAL(toggled(bool)), ag, SLOT(setVisible(bool)));
  ag = new QActionGroup(this);
  ag->addAction(ui->actionMakePath);
  ag->addAction(ui->actionMakeClosedPath);
  ag->addAction(ui->actionDeletePath);
  ag->addAction(ui->actionPathFill);
  ag->setVisible(false);
  connect(ui->actionPath, SIGNAL(toggled(bool)), ag, SLOT(setVisible(bool)));
}

PanelSurface::~PanelSurface()
{
  delete ui;
}

void PanelSurface::ConnectLayer( Layer* layer_in )
{
  PanelLayer::ConnectLayer( layer_in );

  LayerSurface* layer = qobject_cast<LayerSurface*>(layer_in);
  if ( !layer )
  {
    return;
  }

  LayerPropertySurface* p = layer->GetProperty();
  connect( p, SIGNAL(PropertyChanged()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( layer, SIGNAL(SurfaceAnnotationAdded(SurfaceAnnotation*)),
           this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( layer, SIGNAL(SurfaceLabelAdded(SurfaceLabel*)),
           this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( layer, SIGNAL(SurfaceRGBAdded()),
           this, SLOT(UpdateWidgets()), Qt::UniqueConnection );
  connect( layer, SIGNAL(SurfaceLabelDeleted(SurfaceLabel*)), this, SLOT(UpdateWidgets()));
  connect( layer, SIGNAL(SurfaceCurvatureLoaded()), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(SurfaceVectorLoaded()), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(SurfaceOverlayAdded(SurfaceOverlay*)), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(SurfaceSplineAdded(SurfaceSpline*)), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
  connect( layer, SIGNAL(SurfaceSplineDeleted(SurfaceSpline*)), this, SLOT(UpdateWidgets()));

  connect( ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), p, SLOT(SetOpacity(double)) );
  connect( ui->checkBoxShowInfo, SIGNAL(toggled(bool)), p, SLOT(SetShowInfo(bool)) );
  connect( ui->checkBoxShowVertices, SIGNAL(toggled(bool)), p, SLOT(ShowVertices(bool)) );
  connect( ui->colorpickerEdgeColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetEdgeColor(QColor)) );
  connect( ui->colorpickerVectorColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetVectorColor(QColor)) );
  connect( ui->colorpickerSurfaceColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetBinaryColor(QColor)) );
  connect( ui->colorpickerVertexColor, SIGNAL(colorChanged(QColor)), p, SLOT(SetVertexColor(QColor)) );
  connect( ui->comboBoxMeshColor, SIGNAL(currentIndexChanged(int)), p, SLOT(SetMeshColorMap(int)) );
  connect( ui->comboBoxRender, SIGNAL(currentIndexChanged(int)), p, SLOT(SetSurfaceRenderMode(int)) );
  connect( ui->spinBoxEdgeThickness, SIGNAL(valueChanged(int)), p, SLOT(SetEdgeThickness(int)) );
  connect( ui->spinBoxVectorPointSize, SIGNAL(valueChanged(int)), p, SLOT(SetVectorPointSize(int)) );
  connect( ui->spinBoxVertexPointSize, SIGNAL(valueChanged(int)), p, SLOT(SetVertexPointSize(int)) );
  connect( ui->checkBoxUseSurfaceColor2D, SIGNAL(toggled(bool)), p, SLOT(SetUseSurfaceColorOn2D(bool)));
  connect( ui->checkBoxAnnotationOutline, SIGNAL(toggled(bool)), layer, SLOT(SetActiveAnnotationOutline(bool)));
  connect( ui->checkBoxHideIn3DView, SIGNAL(toggled(bool)), layer, SLOT(SetHideIn3D(bool)));
  connect( ui->colorpickerSplineColor, SIGNAL(colorChanged(QColor)), this, SLOT(OnColorPickerSplineColor(QColor)));
  connect(ui->checkBoxSplineProjection, SIGNAL(toggled(bool)), this, SLOT(OnCheckBoxSplineProjection(bool)));
  connect(ui->spinBoxZOrderAnnotation, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxZOrder(int)));
  connect(ui->spinBoxZOrderLabel, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxZOrder(int)));
  connect(ui->spinBoxZOrderOverlay, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxZOrder(int)));

  for (int i = 0; i < layer->GetNumberOfSplines(); i++)
  {
    SurfaceSpline* spline = layer->GetSpline(i);
    connect(spline, SIGNAL(SplineChanged()), this, SLOT(UpdateWidgets()));
  }

  connect(layer, SIGNAL(RGBMapChanged()), this, SLOT(UpdateWidgets()));
  connect(ui->lineEditMappingSurface, SIGNAL(textChanged(QString)), layer, SLOT(SetMappingSurfaceName(QString)));
}

void PanelSurface::DisconnectAllLayers()
{
  PanelLayer::DisconnectAllLayers();
  //    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection("Surface");
  //    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  //    {
  //      LayerSurface* layer = (LayerSurface*)lc->GetLayer(i);
  //      for ( int j = 0; j < allWidgets.size(); j++ )
  //      {
  //        for (int k = 0; k < layer->GetNumberOfSplines(); k++)
  //        {
  //            layer->GetSpline(k)->disconnect( this );
  //            allWidgets[j]->disconnect( layer->GetSpline(k) );
  //        }
  //      }
  //    }
}

void PanelSurface::DoIdle()
{
  // update action status
  BlockAllSignals( true );
  LayerSurface* layer = GetCurrentLayer<LayerSurface*>();
  FSSurface* surf = ( layer ? layer->GetSourceSurface() : NULL );
  ui->actionLockLayer->setEnabled( layer );
  ui->actionLockLayer->setChecked( layer && layer->IsLocked() );
  ui->actionSurfaceMain->setEnabled( layer );
  ui->actionSurfaceMain->setChecked( layer && layer->GetActiveSurface() == FSSurface::SurfaceMain );
  ui->actionSurfaceInflated->setEnabled( surf && surf->IsSurfaceLoaded( FSSurface::SurfaceInflated ) );
  ui->actionSurfaceInflated->setChecked( layer && layer->GetActiveSurface() == FSSurface::SurfaceInflated );
  ui->actionSurfaceOriginal->setEnabled( surf && surf->IsSurfaceLoaded( FSSurface::SurfaceOriginal ) );
  ui->actionSurfaceOriginal->setChecked( layer && layer->GetActiveSurface() == FSSurface::SurfaceOriginal );
  ui->actionSurfaceWhite->setEnabled( surf && surf->IsSurfaceLoaded( FSSurface::SurfaceWhite ) );
  ui->actionSurfaceWhite->setChecked( layer && layer->GetActiveSurface() == FSSurface::SurfaceWhite );
  ui->actionSurfacePial->setEnabled( surf && surf->IsSurfaceLoaded( FSSurface::SurfacePial ) );
  ui->actionSurfacePial->setChecked( layer && layer->GetActiveSurface() == FSSurface::SurfacePial );
  ui->actionShowOverlay->setEnabled(layer && layer->GetNumberOfOverlays() > 0);
  ui->actionShowOverlay->setChecked(layer && layer->GetProperty()->GetShowOverlay());
  ui->actionShowAnnotation->setEnabled(layer && layer->GetNumberOfAnnotations() > 0);
  ui->actionShowAnnotation->setChecked(layer && layer->GetProperty()->GetShowAnnotation());
  ui->actionMoveLayerUp->setEnabled(layer && m_layerCollection
                                    && m_layerCollection->GetLayerIndex(layer) > 0);
  ui->actionMoveLayerDown->setEnabled(layer && m_layerCollection
                                      && m_layerCollection->GetLayerIndex(layer) < m_layerCollection->GetNumberOfLayers()-1);

  int nMode = MainWindow::GetMainWindow()->GetMode();
  ui->actionCut->setChecked(nMode == RenderView::IM_SurfaceCut);
  ui->actionPath->setChecked(nMode == RenderView::IM_SurfacePath);

  ui->actionUndoCut->setEnabled(layer && layer->HasUndoableCut());
  if (nMode != RenderView::IM_SurfacePath)
    m_dlgCustomFill->hide();

  QList<QAction*> acts = m_actGroupSurface->actions();
  foreach (QAction* act, acts)
  {
    if (act->isChecked())
    {
      m_toolButtonSurface->setIcon(act->icon());
      m_toolButtonSurface->setToolTip(act->text());
    }
  }
  BlockAllSignals( false );
}

void PanelSurface::DoUpdateWidgets()
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

  LayerSurface* layer = GetCurrentLayer<LayerSurface*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar &&
         allWidgets[i] != ui->toolbarPath && allWidgets[i]->parentWidget() != ui->toolbarPath)
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  ui->actionCut->setChecked(MainWindow::GetMainWindow()->GetMode() == RenderView::IM_SurfaceCut);
  ui->actionPath->setChecked(MainWindow::GetMainWindow()->GetMode() == RenderView::IM_SurfacePath);
  FSSurface* surf = NULL;
  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->checkBoxHideIn3DView->setChecked(!layer->GetVisibleIn3D());

    surf = layer->GetSourceSurface();
    //    ui->toolbarSurfaces->setVisible(surf->IsSurfaceLoaded( FSSurface::SurfaceOriginal ) || surf->IsSurfaceLoaded( FSSurface::SurfaceInflated ) ||
    //                             surf->IsSurfaceLoaded( FSSurface::SurfaceWhite ) || surf->IsSurfaceLoaded( FSSurface::SurfacePial ) );
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );

    double* rgb = layer->GetProperty()->GetBinaryColor();
    ui->colorpickerSurfaceColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetEdgeColor();
    ui->colorpickerEdgeColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetVectorColor();
    ui->colorpickerEdgeColor->setEnabled(!layer->GetProperty()->GetUseSurfaceColorOn2D());
    ui->checkBoxUseSurfaceColor2D->setChecked(layer->GetProperty()->GetUseSurfaceColorOn2D());
    ui->colorpickerVectorColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ui->lineEditFileName->setText( MyUtils::Win32PathProof(layer->GetFileName()) );
    ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
    ui->spinBoxEdgeThickness->setValue( layer->GetProperty()->GetEdgeThickness() );
    ui->spinBoxVectorPointSize->setValue( layer->GetProperty()->GetVectorPointSize() );

    ui->comboBoxCurvature->setCurrentIndex( layer->GetProperty()->GetCurvatureMap() );
    ChangeLineEditNumber( ui->lineEditMidPoint, layer->GetProperty()->GetThresholdMidPoint() );
    ChangeLineEditNumber( ui->lineEditSlope, layer->GetProperty()->GetThresholdSlope() );
    double range[2];
    layer->GetCurvatureRange( range );
    ui->sliderMidPoint->setValue( (int) ( ( layer->GetProperty()->GetThresholdMidPoint() - range[0] ) / ( range[1] - range[0] ) * 100 ) );
    ui->sliderSlope->setValue( (int) ( layer->GetProperty()->GetThresholdSlope() ) );

    ui->comboBoxRender->setCurrentIndex( layer->GetProperty()->GetSurfaceRenderMode() );
    ui->comboBoxMeshColor->setCurrentIndex( layer->GetProperty()->GetMeshColorMap() );
    ui->checkBoxShowVertices->setChecked( layer->GetProperty()->GetShowVertices() );
    rgb = layer->GetProperty()->GetVertexColor();
    ui->colorpickerVertexColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ChangeSpinBoxValue( ui->spinBoxVertexPointSize, layer->GetProperty()->GetVertexPointSize() );

    double* dPos = layer->GetProperty()->GetPosition();
    ChangeLineEditText( ui->lineEditPositionOffset, QString("%1 %2 %3").arg(dPos[0]).arg(dPos[1]).arg(dPos[2]) );
    ui->checkBoxShowInfo->setChecked( layer->GetProperty()->GetShowInfo() );

    ui->colorpickerSurfaceColor->setVisible(layer->GetActiveRGBMap() < 0);
    QStringList rgb_names = layer->GetRGBMapNames();
    ui->comboBoxColor->clear();
    ui->comboBoxColor->addItem("Solid Color");
    for (int i = 0; i < rgb_names.size(); i++)
      ui->comboBoxColor->addItem(rgb_names[i]);
    ui->comboBoxColor->addItem("Load RGB map...");
    ui->comboBoxColor->setCurrentIndex(layer->GetActiveRGBMap()+1);
  }

  // update vector controls
  ui->comboBoxVectorDisplay->clear();
  ui->comboBoxVectorDisplay->addItem( "Off" );
  if ( surf )
  {
    for ( int i = 0; i < surf->GetNumberOfVectorSets(); i++ )
    {
      ui->comboBoxVectorDisplay->addItem( surf->GetVectorSetName( i ) );
    }
  }
  ui->comboBoxVectorDisplay->addItem( "Load vector data..." );
  ui->comboBoxVectorDisplay->setCurrentIndex( surf ? 1 + surf->GetActiveVector() : 0 );

  // update spline contorls
  //  ui->comboBoxSplineDisplay->clear();
  //  ui->comboBoxSplineDisplay->addItem("Off");
  //  SurfaceSpline* spline = (layer? layer->GetSpline() : NULL);
  //  if (spline && spline->IsValid())
  //  {
  //      ui->comboBoxSplineDisplay->addItem(spline->GetName());
  //  }
  //  ui->comboBoxSplineDisplay->addItem("Load spline data...");
  //  ui->comboBoxSplineDisplay->setCurrentIndex((spline && spline->IsValid() && spline->IsVisible())?1:0);
  //  if (spline)
  //  {
  //    ui->colorpickerSplineColor->setCurrentColor(spline->GetColor());
  //    ui->checkBoxSplineProjection->setChecked(spline->GetProjection());
  //  }

  // update overlay controls
  ui->comboBoxOverlay->clear();
  ui->comboBoxOverlay->addItem( "Off" );
  if ( layer )
  {
    for ( int i = 0; i < layer->GetNumberOfOverlays(); i++ )
    {
      ui->comboBoxOverlay->addItem( layer->GetOverlay( i )->GetName() );
    }
  }
  ui->comboBoxOverlay->addItem( "Load generic..." );
  ui->comboBoxOverlay->addItem( "Load correlation..." );
  ui->comboBoxOverlay->setCurrentIndex( layer ? 1 + layer->GetActiveOverlayIndex() : 0 );
  ShowWidgets(m_widgetsOverlay, layer && layer->GetActiveOverlayIndex() >= 0);
  if ( ui->comboBoxOverlay->currentIndex() == 0 )
  {
    this->m_wndConfigureOverlay->hide();
  }

  // update annotation controls
  ui->comboBoxAnnotation->clear();
  ui->comboBoxAnnotation->addItem( "Off" );
  ui->comboBoxAnnotation->setItemData(0, "Off", Qt::ToolTipRole);
  if ( layer )
  {
    for ( int i = 0; i < layer->GetNumberOfAnnotations(); i++ )
    {
      ui->comboBoxAnnotation->addItem( layer->GetAnnotation( i )->GetName() );
      ui->comboBoxAnnotation->setItemData(i+1, layer->GetAnnotation(i)->GetFilename(), Qt::ToolTipRole);
    }
  }
  ui->comboBoxAnnotation->addItem( "Load from file..." );
  ui->comboBoxAnnotation->setItemData(ui->comboBoxAnnotation->count()-1, "Load from file", Qt::ToolTipRole);
  ui->comboBoxAnnotation->setCurrentIndex( layer ? 1 + layer->GetActiveAnnotationIndex() : 0 );

  // update label controls
  QList<SurfaceLabel*> selected_labels;
  if (layer && ui->treeWidgetLabels->topLevelItemCount() == layer->GetNumberOfLabels())
    selected_labels = GetSelectedLabels();
  ui->treeWidgetLabels->clear();
  if (layer)
  {
    for (int i = 0; i < layer->GetNumberOfLabels(); i++)
    {
      SurfaceLabel* label = layer->GetLabel(i);
      QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetLabels);
      item->setText(0, label->GetName());
      item->setData(0, Qt::UserRole, QVariant::fromValue((QObject*)label));
      item->setCheckState(0, label->IsVisible() ? Qt::Checked : Qt::Unchecked);
      item->setFlags(item->flags() | Qt::ItemIsEditable);
      if (layer->GetActiveLabel() == label)
        ui->treeWidgetLabels->setCurrentItem(item);
      if (selected_labels.contains(label))
      {
        item->setSelected(true);
      }
    }
    ui->pushButtonDeleteLabel->setEnabled(layer->GetNumberOfLabels() > 0);
  }

  // update spline contorls
  QList<SurfaceSpline*> selected_splines;
  if (layer && ui->treeWidgetSplines->topLevelItemCount() == layer->GetNumberOfSplines())
    selected_splines = GetSelectedSplines();
  ui->treeWidgetSplines->clear();
  if (layer)
  {
    for (int i = 0; i < layer->GetNumberOfSplines(); i++)
    {
      SurfaceSpline* spline = layer->GetSpline(i);
      QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetSplines);
      item->setText(0, spline->GetName());
      item->setData(0, Qt::UserRole, QVariant::fromValue(reinterpret_cast<quintptr>(spline)));
      item->setCheckState(0, spline->IsVisible() ? Qt::Checked : Qt::Unchecked);
      if (layer->GetActiveSpline() == spline)
        ui->treeWidgetSplines->setCurrentItem(item);
      if (selected_splines.contains(spline))
      {
        item->setSelected(true);
      }
    }
    ui->pushButtonDeleteSpline->setEnabled(layer->GetNumberOfSplines() > 0);
  }

  if ( layer && layer->GetActiveAnnotation() )
  {
    ui->checkBoxAnnotationOutline->setChecked(layer->GetActiveAnnotation()->GetShowOutline());
  }

  int nCurvatureMap = layer ? layer->GetProperty()->GetCurvatureMap() : 0;
  ShowWidgets( m_widgetsMidPoint, nCurvatureMap != LayerPropertySurface::CM_Off );
  ShowWidgets( m_widgetsSlope,    nCurvatureMap == LayerPropertySurface::CM_Threshold );
  ShowWidgets( m_widgetsVector,   ui->comboBoxVectorDisplay->currentIndex() > 0 );
  ShowWidgets( m_widgetsVertex,   ui->checkBoxShowVertices->isChecked() );
  ShowWidgets( m_widgetsMesh,     layer && layer->GetProperty()->GetSurfaceRenderMode() != LayerPropertySurface::SM_Surface );
  ShowWidgets( m_widgetsLabel,    layer && layer->GetActiveLabelIndex() >= 0 );
  ShowWidgets( m_widgetsSpline,   layer && layer->GetActiveSpline());
  ShowWidgets( m_widgetsAnnotation, layer && layer->GetActiveAnnotation());
  UpdateLabelWidgets(false);
  ui->colorpickerSurfaceColor->setEnabled( layer ); // && nCurvatureMap != LayerPropertySurface::CM_Threshold );

  bool isInflated = (layer && layer->GetFileName().contains("inflated"));
  ui->labelMapCursorTo->setVisible(isInflated);
  ui->lineEditMappingSurface->setVisible(isInflated);
  if (isInflated)
    ui->lineEditMappingSurface->setText(layer->GetMappingSurfaceName());

  if (layer)
  {
    ui->spinBoxZOrderAnnotation->setValue(layer->GetProperty()->GetZOrderAnnotation());
    ui->spinBoxZOrderLabel->setValue(layer->GetProperty()->GetZOrderLabel());
    ui->spinBoxZOrderOverlay->setValue(layer->GetProperty()->GetZOrderOverlay());
  }

  UpdateSplineWidgets();

  BlockAllSignals( false );
}

void PanelSurface::UpdateLabelWidgets(bool block_signals)
{
  LayerSurface* layer = GetCurrentLayer<LayerSurface*>();
  if (layer && layer->GetActiveLabel())
  {
    if (block_signals)
      BlockAllSignals(true);

    SurfaceLabel* label = layer->GetActiveLabel();
    double* rgb = label->GetColor();
    ui->colorpickerLabelColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ui->checkBoxLabelOutline->setChecked(label->GetShowOutline());

    ChangeLineEditNumber(ui->lineEditLabelThreshold, label->GetThreshold());
    ChangeLineEditNumber(ui->lineEditLabelHeatscaleMin, label->GetHeatscaleMin());
    ChangeLineEditNumber(ui->lineEditLabelHeatscaleMax, label->GetHeatscaleMax());
    ChangeLineEditNumber(ui->lineEditLabelOpacity, label->GetOpacity());
    ui->comboBoxLabelColorCode->setCurrentIndex(label->GetColorCode());
    ui->labelHeatscaleRange->setVisible(label->GetColorCode() == SurfaceLabel::Heatscale);
    ui->lineEditLabelHeatscaleMin->setVisible(label->GetColorCode() == SurfaceLabel::Heatscale);
    ui->lineEditLabelHeatscaleMax->setVisible(label->GetColorCode() == SurfaceLabel::Heatscale);
    ui->colorpickerLabelColor->setVisible(label->GetColorCode() == SurfaceLabel::SolidColor);

    ui->toolButtonLabelUp->setEnabled(label != layer->GetLabel(0));
    ui->toolButtonLabelDown->setEnabled(label != layer->GetLabel(layer->GetNumberOfLabels()-1));

    QMenu* menu = new QMenu();
    ui->toolButtonMoreLabelOptions->setMenu(menu);
    QAction* act = new QAction("Go To Centroid", this);
    connect(act, SIGNAL(triggered()), MainWindow::GetMainWindow(), SLOT(OnGoToSurfaceLabel()));
    menu->addAction(act);
    act = new QAction("Resample", this);
    connect(act, SIGNAL(triggered()), this, SLOT(OnLabelResample()));
    menu->addAction(act);
    act = new QAction("Dilate/Erode/Open/Close...", this);
    connect(act, SIGNAL(triggered()), this, SLOT(OnLabelMoreOps()));
    menu->addAction(act);
    menu->addSeparator();
    act = new QAction("Save As...", this);
    connect(act, SIGNAL(triggered()), this, SLOT(OnSaveLabelAs()));
    menu->addAction(act);

    if (block_signals)
      BlockAllSignals(false);
  }
}

void PanelSurface::UpdateSplineWidgets()
{
  LayerSurface* layer = GetCurrentLayer<LayerSurface*>();
  if (layer && layer->GetActiveSpline())
  {
    BlockAllSignals(true);

    SurfaceSpline* spline = layer->GetActiveSpline();
    ui->colorpickerSplineColor->setCurrentColor( spline->GetColor() );
    ui->checkBoxSplineProjection->setChecked(spline->GetProjection());

    BlockAllSignals(false);
  }
}

QList<SurfaceLabel*> PanelSurface::GetSelectedLabels()
{
  QList<SurfaceLabel*> selected_labels;
  QList<QTreeWidgetItem*> selected_items = ui->treeWidgetLabels->selectedItems();
  foreach (QTreeWidgetItem* item, selected_items)
  {
    selected_labels << qobject_cast<SurfaceLabel*>(item->data(0, Qt::UserRole).value<QObject*>());
  }
  return selected_labels;
}

QList<SurfaceSpline*> PanelSurface::GetSelectedSplines()
{
  QList<SurfaceSpline*> selected_splines;
  QList<QTreeWidgetItem*> selected_items = ui->treeWidgetSplines->selectedItems();
  foreach (QTreeWidgetItem* item, selected_items)
  {
    selected_splines << reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  }
  return selected_splines;
}


void PanelSurface::OnChangeSurfaceType( QAction* act )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    layer->SetActiveSurface( act->data().toInt() );
  }
}

void PanelSurface::OnSliderOpacity( int nVal )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    layer->GetProperty()->SetOpacity( nVal / 100.0 );
  }
}

void PanelSurface::OnSliderMidPoint( int nVal )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    double range[2];
    layer->GetCurvatureRange( range );
    layer->GetProperty()->SetThresholdMidPoint( nVal / 100.0 * ( range[1] - range[0] ) + range[0] );
  }
}

void PanelSurface::OnSliderSlope( int nVal )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    double fMin = 0;
    double fMax = 100;
    layer->GetProperty()->SetThresholdSlope( nVal / 100.0 * ( fMax - fMin ) + fMin );
  }
}

void PanelSurface::OnComboCurvature( int nSel )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    if ( nSel < 3 )
    {
      layer->GetProperty()->SetCurvatureMap( nSel );
    }
    else
    {
      // load new curvature maps
      MainWindow::GetMainWindow()->LoadSurfaceCurvature();
    }
  }
}

void PanelSurface::OnLineEditMidPoint( const QString& text )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    bool bOK;
    double dval = text.toDouble( &bOK );
    if ( layer && dval && dval != layer->GetProperty()->GetThresholdMidPoint() )
    {
      layer->GetProperty()->SetThresholdMidPoint( dval );
    }
  }
}

void PanelSurface::OnLineEditSlope( const QString& text )
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    bool bOK;
    double dval = text.toDouble( &bOK );
    if ( layer && dval && dval != layer->GetProperty()->GetThresholdSlope() )
    {
      layer->GetProperty()->SetThresholdSlope( dval );
    }
  }
}

void PanelSurface::OnLineEditLabelThreshold(const QString &text)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  bool bOK;
  double dval = text.toDouble( &bOK );
  /*
  if (surf && surf->GetActiveLabel() && bOK)
  {
    surf->GetActiveLabel()->SetThreshold(dval);
  }
  */
  if (surf && bOK)
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetThreshold(dval);
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnComboOverlay( int nSel_in )
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    int nSel = nSel_in - 1;
    if ( nSel < surf->GetNumberOfOverlays() )
    {
      surf->SetActiveOverlay( nSel );
    }
    else
    {
      // load new overlay map
      MainWindow::GetMainWindow()->LoadSurfaceOverlay(nSel > surf->GetNumberOfOverlays());
    }
    UpdateWidgets();
  }
}

void PanelSurface::OnComboAnnotation( int nSel_in )
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    int nSel = nSel_in - 1;
    if ( nSel < surf->GetNumberOfAnnotations() )
    {
      surf->SetActiveAnnotation( nSel );
    }
    else
    {
      // load new overlay map
      MainWindow::GetMainWindow()->LoadSurfaceAnnotation();
    }
    UpdateWidgets();
  }
}

void PanelSurface::OnButtonLoadLabel()
{
  MainWindow::GetMainWindow()->LoadSurfaceLabel();
  UpdateWidgets();
}

void PanelSurface::OnButtonDeleteLabel()
{
  QTreeWidgetItem* item = ui->treeWidgetLabels->currentItem();
  if (!item)
    return;

  SurfaceLabel* label = qobject_cast<SurfaceLabel*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( label )
  {
    LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
    if ( surf )
    {
      surf->DeleteLabel(label);
    }
    UpdateLabelWidgets();
  }
}

void PanelSurface::OnButtonNewLabel()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    surf->CreateNewLabel();
  }
}

void PanelSurface::OnButtonSaveLabel()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if (surf)
  {
    SurfaceLabel* label = surf->GetActiveLabel();
    if (label)
    {
      QDir dir = QFileInfo(surf->GetFileName()).absoluteDir();
      dir.cd("label");
      QString fn = label->GetFileName();
      if (fn.isEmpty())
      {
        QString def_fn = label->GetName() + ".label";
        def_fn = dir.absoluteFilePath(def_fn);
        fn = QFileDialog::getSaveFileName( this, "Select label file",
                                           def_fn,
                                           "Label files (*)");
      }
      if (!fn.isEmpty())
      {
        label->SaveToFile(fn);
      }
    }
  }
}

void PanelSurface::OnSaveLabelAs()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if (surf)
  {
    SurfaceLabel* label = surf->GetActiveLabel();
    if (label)
    {
      QString fn = label->GetFileName();
      if (fn.isEmpty())
      {
        OnButtonSaveLabel();
      }
      else
      {
        QString def_fn = label->GetName() + ".label";
        def_fn = QFileInfo(fn).absolutePath() + "/" + def_fn;
        fn = QFileDialog::getSaveFileName( this, "Select label file",
                                           def_fn,
                                           "Label files (*)");
        if (!fn.isEmpty())
        {
          label->SaveToFile(fn);
        }
      }
    }
  }
}

void PanelSurface::OnComboVector( int nSel )
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    int nVector = nSel - 1;
    if ( nVector < surf->GetNumberOfVectorSets() )
    {
      surf->SetActiveVector( nVector );
    }
    else
    {
      // load new
      MainWindow::GetMainWindow()->LoadSurfaceVector();
    }
    UpdateWidgets();
  }
}

void PanelSurface::OnComboSpline(int nSel)
{
  Q_UNUSED(nSel);
  //  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  //  if ( surf )
  //  {
  //    if (surf->GetSpline()->IsValid() && nSel == 1)
  //      surf->GetSpline()->SetVisible(true);
  //    else if (nSel == 0)
  //      surf->GetSpline()->SetVisible(false);
  //    else
  //    {
  //      MainWindow::GetMainWindow()->LoadSurfaceSpline();
  //    }
  //    UpdateWidgets();
  //  }
}

void PanelSurface::OnButtonConfigureOverlay()
{
  m_wndConfigureOverlay->show();
}

void PanelSurface::OnButtonRemoveOverlay()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    surf->RemoveCurrentOverlay();
  }
  UpdateWidgets();
}

void PanelSurface::OnEditPositionOffset()
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    QStringList args = ui->lineEditPositionOffset->text().trimmed().split(" ", QString::SkipEmptyParts);
    args << "n/a" << "n/a" << "n/a";
    bool bOK;
    double pos[3];
    pos[0] = args[0].toDouble(&bOK);
    if (bOK)
    {
      pos[1] = args[1].toDouble(&bOK);
    }
    if (bOK)
    {
      pos[2] = args[2].toDouble(&bOK);
    }
    if (bOK)
    {
      layer->GetProperty()->SetPosition( pos );
    }
    else
    {
      //  QMessageBox::information(this, "Error", "Please enter 3 values for position offset.");
    }
  }
}

void PanelSurface::OnLabelItemChanged(QTreeWidgetItem *item)
{
  SurfaceLabel* label = qobject_cast<SurfaceLabel*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( label )
  {
    label->SetVisible( item->checkState( 0 ) == Qt::Checked );
    if (item->text(0) != label->GetName())
      label->SetName(item->text(0));
  }
}

void PanelSurface::OnLabelItemDoubleClicked(QTreeWidgetItem *item)
{
  SurfaceLabel* label = qobject_cast<SurfaceLabel*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( label )
  {
    LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
    if ( surf )
    {
      surf->MoveLabelToTop(label);
      UpdateWidgets();
    }
  }
}

void PanelSurface::OnCurrentLabelItemChanged(QTreeWidgetItem *item)
{
  SurfaceLabel* label = qobject_cast<SurfaceLabel*>(item->data( 0, Qt::UserRole ).value<QObject*>());
  if ( label )
  {
    LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
    if ( surf )
    {
      surf->blockSignals(true);
      surf->SetActiveLabel(label);
      surf->blockSignals(false);
      UpdateLabelWidgets();
    }
  }
}

void PanelSurface::OnSplineItemChanged(QTreeWidgetItem *item)
{
  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  if ( spline )
  {
    spline->SetVisible( item->checkState( 0 ) == Qt::Checked );
  }
}

void PanelSurface::OnCurrentSplineItemChanged(QTreeWidgetItem *item)
{
  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  if ( spline )
  {
    LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
    if ( surf )
    {
      surf->blockSignals(true);
      surf->SetActiveSpline(spline);
      surf->blockSignals(false);
      UpdateSplineWidgets();
    }
  }
}

void PanelSurface::OnToggleOverlay(bool bShow)
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    layer->GetProperty()->SetShowOverlay(bShow);
  }
}

void PanelSurface::OnToggleAnnotation(bool bShow)
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    layer->GetProperty()->SetShowAnnotation(bShow);
  }
}

void PanelSurface::OnComboLabelColorCode(int nSel)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->SetColorCode(nSel);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnLineEditLabelHeatscaleMin(const QString &text)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  bool bOK;
  double dval = text.toDouble( &bOK );
  if (surf && bOK)
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetHeatscaleMin(dval);
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnLineEditLabelHeatscaleMax(const QString &text)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  bool bOK;
  double dval = text.toDouble( &bOK );
  if (surf && bOK)
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetHeatscaleMax(dval);
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnCheckBoxLabelOutline(bool outline)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetShowOutline(outline);
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnColorPickerLabelColor(const QColor &color)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetColor(color.redF(), color.greenF(), color.blueF());
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnLockLayer(bool b)
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    layer->Lock(b);
  }
}

void PanelSurface::OnComboColor(int nSel)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    if (nSel == 0)
      surf->SetActiveRGBMap(-1);
    else if (nSel > surf->GetNumberOfRGBMaps())
    {
      if (!MainWindow::GetMainWindow()->LoadSurfaceRGBMap())
        UpdateWidgets();
    }
    else
      surf->SetActiveRGBMap(nSel-1);
  }
}

void PanelSurface::OnButtonLoadSpline()
{
  MainWindow::GetMainWindow()->LoadSurfaceSpline();
  UpdateWidgets();
}

void PanelSurface::OnButtonDeleteSpline()
{
  QTreeWidgetItem* item = ui->treeWidgetSplines->currentItem();
  if (!item)
    return;

  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  if ( spline )
  {
    LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
    if ( surf )
    {
      surf->DeleteSpline(spline);
    }
  }
}

void PanelSurface::OnColorPickerSplineColor(const QColor &color)
{
  QTreeWidgetItem* item = ui->treeWidgetSplines->currentItem();
  if (!item)
    return;

  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  if ( spline )
  {
    spline->SetColor(color);
  }
}

void PanelSurface::OnCheckBoxSplineProjection(bool b)
{
  QTreeWidgetItem* item = ui->treeWidgetSplines->currentItem();
  if (!item)
    return;

  SurfaceSpline* spline = reinterpret_cast<SurfaceSpline*>(item->data(0, Qt::UserRole).value<quintptr>());
  if ( spline )
  {
    spline->SetProjection(b);
  }
}

void PanelSurface::OnButtonEditCut(bool b)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (b)
  {
    mainwnd->SetMode(RenderView::IM_SurfaceCut);
    ui->actionPath->setChecked(false);
  }
  else if (mainwnd->GetMode() != RenderView::IM_SurfacePath)
    mainwnd->SetMode(RenderView::IM_Navigate);
}

void PanelSurface::OnButtonEditPath(bool b)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if (b)
  {
    mainwnd->SetMode(RenderView::IM_SurfacePath);
    ui->actionCut->setChecked(false);
  }
  else if (mainwnd->GetMode() != RenderView::IM_SurfaceCut)
    mainwnd->SetMode(RenderView::IM_Navigate);
}

void PanelSurface::OnButtonCutLine()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetMarks())
  {
    if (surf->GetMarks()->GetPathVerts().size() < 2)
      QMessageBox::warning(this->parentWidget(), "Error", "Need at least 2 marks to cut line");
    else
      surf->GetMarks()->MakeCutLine(false);
  }
}

void PanelSurface::OnButtonCutClosedLine()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetMarks())
  {
    if (surf->GetMarks()->GetPathVerts().size() < 3)
      QMessageBox::warning(this->parentWidget(), "Error", "Need at least 3 marks to cut closed line");
    else
      surf->GetMarks()->MakeCutLine(true);
  }
}

void PanelSurface::OnButtonClearCuts()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf)
  {
    surf->ClearAllCuts();
  }
}

void PanelSurface::OnButtonFillUncutArea()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf)
  {
    int vno = surf->GetCurrentVertex();
    if (vno < 0)
      QMessageBox::information(this->parentWidget(), "Fill Cut Area", "Please move cursor to a valid vertex");
    else if (surf->IsVertexRipped(vno))
      QMessageBox::information(this->parentWidget(), "Fill Cut Area", "Please move cursor to an uncut vertex");
    else
    {
      surf->FillUncutArea(vno);
      //  ui->actionCut->setChecked(false);
    }
  }
}

void PanelSurface::OnButtonUndoCut()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf)
    surf->UndoCut();
}

void PanelSurface::OnLabelResample()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->GetActiveLayer("MRI"));
  if (!mri)
  {
    QMessageBox::warning(this, "Error", "Could not find a MRI template for resampling");
    return;
  }
  if ( surf && surf->GetActiveLabel())
    surf->GetActiveLabel()->Resample(mri);
}

void PanelSurface::OnLabelMoreOps()
{
  m_dlgLabelOps->show();
  m_dlgLabelOps->raise();
}

void PanelSurface::OnLabelOperation(const QVariantMap &op)
{
  QString op_str = op["operation"].toString();
  int nTimes = op["times"].toInt();
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetActiveLabel())
  {
    SurfaceLabel* label = surf->GetActiveLabel();
    if (op_str == "dilate")
      label->Dilate(nTimes);
    else if (op_str == "erode")
      label->Erode(nTimes);
    else if (op_str == "open")
      label->Open(nTimes);
    else if (op_str == "close")
      label->Close(nTimes);
  }
}

void PanelSurface::OnSpinBoxZOrder(int nOrder)
{
  QList<LayerSurface*> layers = GetSelectedLayers<LayerSurface*>();
  foreach (LayerSurface* layer, layers)
  {
    if (sender() == ui->spinBoxZOrderAnnotation)
      layer->GetProperty()->SetZOrderAnnotation(nOrder);
    else if (sender() == ui->spinBoxZOrderLabel)
      layer->GetProperty()->SetZOrderLabel(nOrder);
    else if (sender() == ui->spinBoxZOrderOverlay)
      layer->GetProperty()->SetZOrderOverlay(nOrder);
  }
}

void PanelSurface::OnButtonMakePath()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetMarks())
  {
    if (surf->GetMarks()->GetPathVerts().size() < 2)
      QMessageBox::warning(this->parentWidget(), "Error", "Need at least 2 marks to make path");
    else
      surf->GetMarks()->MakePath(false);
  }
}

void PanelSurface::OnButtonMakeClosedPath()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetMarks())
  {
    if (surf->GetMarks()->GetPathVerts().size() < 3)
      QMessageBox::warning(this->parentWidget(), "Error", "Need at least 3 marks to make closed path");
    else
      surf->GetMarks()->MakePath(true);
  }
}

void PanelSurface::OnButtonDeletePath()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf && surf->GetActivePath())
  {
    surf->DeleteActivePath();
  }
}

void PanelSurface::OnButtonCustomFillPath()
{
  m_dlgCustomFill->show();
  m_dlgCustomFill->raise();
}

void PanelSurface::OnCustomFillTriggered(const QVariantMap &options)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf)
  {
    int vno = surf->GetLastMark();
    if (vno < 0)
      vno = surf->GetCurrentVertex();
    if (vno < 0)
      QMessageBox::information(this->parentWidget(), "Fill Cut Area", "Please move cursor to a valid vertex");
    else if (surf->IsVertexRipped(vno))
      QMessageBox::information(this->parentWidget(), "Fill Cut Area", "Please move cursor to an uncut vertex");
    else
      surf->FillPath(vno, options);
  }
}

void PanelSurface::OnButtonClearMarks()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if (surf)
    surf->ClearMarks();
}

void PanelSurface::OnLineEditLabelOpacity(const QString &text)
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  bool bOK;
  double dval = text.toDouble( &bOK );
  if (surf && bOK && dval >= 0 && dval <= 1)
  {
    QList<SurfaceLabel*> labels = GetSelectedLabels();
    foreach (SurfaceLabel* label, labels)
    {
      label->blockSignals(true);
      label->SetOpacity(dval);
      label->blockSignals(false);
    }
    surf->UpdateColorMap();
  }
}

void PanelSurface::OnButtonLabelUp()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  QTreeWidgetItem* curItem = ui->treeWidgetLabels->currentItem();
  if (curItem && surf)
  {
    int n = ui->treeWidgetLabels->indexOfTopLevelItem(curItem);
    if (n > 0)
    {
      ui->treeWidgetLabels->takeTopLevelItem(n);
      ui->treeWidgetLabels->insertTopLevelItem(n-1, curItem);
      surf->MoveLabelUp(surf->GetLabel(n));
      ui->treeWidgetLabels->setCurrentItem(curItem);
    }
  }
}

void PanelSurface::OnButtonLabelDown()
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  QTreeWidgetItem* curItem = ui->treeWidgetLabels->currentItem();
  if (curItem)
  {
    int n = ui->treeWidgetLabels->indexOfTopLevelItem(curItem);
    if (n < ui->treeWidgetLabels->topLevelItemCount()-1)
    {
      ui->treeWidgetLabels->takeTopLevelItem(n);
      ui->treeWidgetLabels->insertTopLevelItem(n+1, curItem);
      surf->MoveLabelDown(surf->GetLabel(n));
      ui->treeWidgetLabels->setCurrentItem(curItem);
    }
  }
}
