/**
 * @file  PanelSurface.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/11/01 19:21:06 $
 *    $Revision: 1.47 $
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

PanelSurface::PanelSurface(QWidget *parent) :
  PanelLayer(parent),
  ui(new Ui::PanelSurface)
{
  ui->setupUi(this);
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  ui->toolbar->insertAction(ui->actionSurfaceMain, mainwnd->ui->actionLoadSurface);
  ui->toolbar->insertAction(ui->actionSurfaceMain, mainwnd->ui->actionCloseSurface);
  ui->toolbar->insertSeparator(ui->actionSurfaceMain);

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
                 << ui->checkBoxLabelOutline;

  m_widgetsSpline << ui->colorpickerSplineColor
                  << ui->labelSplineColor
                  << ui->checkBoxSplineProjection;

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

  LayerCollection* lc = mainwnd->GetLayerCollection("Surface");
  PanelLayer::InitializeLayerList( ui->treeWidgetLayers, lc );
  connect( ui->actionLockLayer, SIGNAL(toggled(bool)), lc, SLOT(LockCurrent(bool)) );

  m_wndConfigureOverlay = new WindowConfigureOverlay( this );
  m_wndConfigureOverlay->hide();
  connect( mainwnd->GetLayerCollection("Surface"), SIGNAL(ActiveLayerChanged(Layer*)),
           m_wndConfigureOverlay, SLOT(OnActiveSurfaceChanged(Layer*)));
  connect(m_wndConfigureOverlay, SIGNAL(ActiveFrameChanged()), mainwnd, SLOT(UpdateInfoPanel()));
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
  connect( layer, SIGNAL(SurfaceCurvatureLoaded()), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(SurfaceVectorLoaded()), this, SLOT(UpdateWidgets()) );
  connect( layer, SIGNAL(SurfaceOverlayAdded(SurfaceOverlay*)), this, SLOT(UpdateWidgets()) );
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
  connect( ui->colorpickerLabelColor, SIGNAL(colorChanged(QColor)), layer, SLOT(SetActiveLabelColor(QColor)));
  connect( ui->checkBoxLabelOutline, SIGNAL(toggled(bool)), layer, SLOT(SetActiveLabelOutline(bool)));
  connect( ui->checkBoxAnnotationOutline, SIGNAL(toggled(bool)), layer, SLOT(SetActiveAnnotationOutline(bool)));

  SurfaceSpline* spline = layer->GetSpline();
  connect( ui->colorpickerSplineColor, SIGNAL(colorChanged(QColor)), spline, SLOT(SetColor(QColor)));
  connect(ui->checkBoxSplineProjection, SIGNAL(toggled(bool)), spline, SLOT(SetProjection(bool)));
  connect(spline, SIGNAL(SplineChanged()), this, SLOT(UpdateWidgets()));
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
  BlockAllSignals( false );
}

void PanelSurface::DoUpdateWidgets()
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

  LayerSurface* layer = GetCurrentLayer<LayerSurface*>();
  for ( int i = 0; i < this->allWidgets.size(); i++ )
  {
    if ( allWidgets[i] != ui->toolbar && allWidgets[i]->parentWidget() != ui->toolbar )
    {
      allWidgets[i]->setEnabled(layer);
    }
  }
  FSSurface* surf = NULL;
  ui->lineEditFileName->clear();
  if ( layer )
  {
    ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
    ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );

    double* rgb = layer->GetProperty()->GetBinaryColor();
    ui->colorpickerSurfaceColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetEdgeColor();
    ui->colorpickerEdgeColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    rgb = layer->GetProperty()->GetVectorColor();
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

    surf = layer->GetSourceSurface();
    ui->comboBoxRender->setCurrentIndex( layer->GetProperty()->GetSurfaceRenderMode() );
    ui->comboBoxMeshColor->setCurrentIndex( layer->GetProperty()->GetMeshColorMap() );
    ui->checkBoxShowVertices->setChecked( layer->GetProperty()->GetShowVertices() );
    rgb = layer->GetProperty()->GetVertexColor();
    ui->colorpickerVertexColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ChangeSpinBoxValue( ui->spinBoxVertexPointSize, layer->GetProperty()->GetVertexPointSize() );

    double* dPos = layer->GetProperty()->GetPosition();
    ChangeLineEditText( ui->lineEditPositionOffset, QString("%1 %2 %3").arg(dPos[0]).arg(dPos[1]).arg(dPos[2]) );
    ui->checkBoxShowInfo->setChecked( layer->GetProperty()->GetShowInfo() );
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
  ui->comboBoxSplineDisplay->clear();
  ui->comboBoxSplineDisplay->addItem("Off");
  SurfaceSpline* spline = (layer? layer->GetSpline() : NULL);
  if (spline && spline->IsValid())
  {
      ui->comboBoxSplineDisplay->addItem(spline->GetName());
  }
  ui->comboBoxSplineDisplay->addItem("Load spline data...");
  ui->comboBoxSplineDisplay->setCurrentIndex((spline && spline->IsValid() && spline->IsVisible())?1:0);
  if (spline)
  {
    ui->colorpickerSplineColor->setCurrentColor(spline->GetColor());
    ui->checkBoxSplineProjection->setChecked(spline->GetProjection());
  }

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
  ui->pushButtonConfigureOverlay->setVisible( layer && layer->GetActiveOverlayIndex() >= 0 );
  if ( ui->comboBoxOverlay->currentIndex() == 0 )
  {
    this->m_wndConfigureOverlay->hide();
  }

  // update annotation controls
  ui->comboBoxAnnotation->clear();
  ui->comboBoxAnnotation->addItem( "Off" );
  if ( layer )
  {
    for ( int i = 0; i < layer->GetNumberOfAnnotations(); i++ )
    {
      ui->comboBoxAnnotation->addItem( layer->GetAnnotation( i )->GetName() );
    }
  }
  ui->comboBoxAnnotation->addItem( "Load from file..." );
  ui->comboBoxAnnotation->setCurrentIndex( layer ? 1 + layer->GetActiveAnnotationIndex() : 0 );

  // update label controls
  ui->comboBoxLabel->clear();
  if ( layer )
  {
    ui->comboBoxLabel->addItem( "Off" );
    if ( layer->GetNumberOfLabels() > 0 )
    {
      for ( int i = 0; i < layer->GetNumberOfLabels(); i++ )
      {
        ui->comboBoxLabel->addItem( layer->GetLabel( i )->GetName() );
      }
    }
  }
  ui->comboBoxLabel->addItem( "Load from file..." );
  if ( layer && layer->GetActiveLabelIndex() >= 0 )
  {
    ui->comboBoxLabel->setCurrentIndex( layer->GetActiveLabelIndex()+1 );
  }
  else
  {
    ui->comboBoxLabel->setCurrentIndex( 0 );
  }
  if ( layer && layer->GetActiveLabel() )
  {
    double* rgb = layer->GetActiveLabel()->GetColor();
    ui->colorpickerLabelColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
    ui->checkBoxLabelOutline->setChecked(layer->GetActiveLabel()->GetShowOutline());
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
  ShowWidgets(m_widgetsSpline, spline && spline->IsValid() && spline->IsVisible());
  ui->checkBoxAnnotationOutline->setVisible(layer && layer->GetActiveAnnotation());
  ui->colorpickerSurfaceColor->setEnabled( layer ); // && nCurvatureMap != LayerPropertySurface::CM_Threshold );

  BlockAllSignals( false );
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

void PanelSurface::OnComboLabel( int nSel_in )
{
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    int nSel = nSel_in - 1;
    if ( nSel >= surf->GetNumberOfLabels() )
    {
      MainWindow::GetMainWindow()->LoadSurfaceLabel();
    }
    else
    {
      surf->SetActiveLabel( nSel );
    }
    UpdateWidgets();
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
  LayerSurface* surf = GetCurrentLayer<LayerSurface*>();
  if ( surf )
  {
    if (surf->GetSpline()->IsValid() && nSel == 1)
      surf->GetSpline()->SetVisible(true);
    else if (nSel == 0)
      surf->GetSpline()->SetVisible(false);
    else
    {
      MainWindow::GetMainWindow()->LoadSurfaceSpline();
    }
    UpdateWidgets();
  }
}

void PanelSurface::OnButtonConfigureOverlay()
{
  m_wndConfigureOverlay->show();
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
      QMessageBox::information(this, "Error", "Please enter 3 values for position offset.");
    }
  }
}
