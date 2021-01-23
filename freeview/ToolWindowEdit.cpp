/*
 * Original Author: Ruopeng Wang
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
#include "ToolWindowEdit.h"
#include "ui_ToolWindowEdit.h"
#include "Interactor2DVoxelEdit.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "Contour2D.h"
#include "BrushProperty.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "RenderView2D.h"
#include "DialogReplaceLabel.h"
#include <QTimer>
#include <QSettings>
#include <QDebug>
#include "LayerPropertyMRI.h"
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif

ToolWindowEdit::ToolWindowEdit(QWidget *parent) :
  QWidget(parent),
  UIUpdateHelper(),
  ui(new Ui::ToolWindowEdit)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );
#ifndef Q_OS_MAC
  ui->line->hide();
#endif
  QActionGroup* ag = new QActionGroup( this );
  ag->addAction( ui->actionContour );
  ag->addAction( ui->actionFill );
  ag->addAction( ui->actionFreeHand );
  ag->addAction( ui->actionLiveWire );
  ag->addAction( ui->actionPolyLine );
  ag->addAction( ui->actionColorPicker );
  ag->addAction( ui->actionClone );
  ag->addAction( ui->actionAutoSeg);
  ag->addAction( ui->actionShift );
  ag->setExclusive( true );

  ui->actionContour->setData( Interactor2DVoxelEdit::EM_Contour );
  ui->actionColorPicker->setData( Interactor2DVoxelEdit::EM_ColorPicker );
  ui->actionFill->setData( Interactor2DVoxelEdit::EM_Fill );
  ui->actionFreeHand->setData( Interactor2DVoxelEdit::EM_Freehand );
  ui->actionLiveWire->setData( Interactor2DVoxelEdit::EM_Livewire );
  ui->actionPolyLine->setData( Interactor2DVoxelEdit::EM_Polyline );
  ui->actionClone->setData( Interactor2DVoxelEdit::EM_Clone );
  ui->actionAutoSeg->setData( Interactor2DVoxelEdit::EM_GeoSeg);
  ui->actionShift->setData( Interactor2DVolumeEdit::EM_Shift );
  ui->colorPickerGeoInside->setCurrentColor(Qt::green);
  ui->colorPickerGeoOutside->setCurrentColor(Qt::red);
  ui->colorPickerGeoFill->setCurrentColor(Qt::yellow);
  ui->widgetBusyIndicator->setColor(Qt::darkGray);
  ui->widgetBusyIndicator->setFixedSize(QSize(20,20));
  connect( ag, SIGNAL(triggered(QAction*)), this, SLOT(OnEditMode(QAction*)) );
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  BrushProperty* bp = mainwnd->GetBrushProperty();
  connect(ui->spinBoxBrushSize, SIGNAL(valueChanged(int)), bp, SLOT(SetBrushSize(int)));
  connect(ui->spinBoxTolerance, SIGNAL(valueChanged(int)), bp, SLOT(SetBrushTolerance(int)));
  connect(ui->comboBoxReference, SIGNAL(currentIndexChanged(int)), this, SLOT(OnComboReference(int)));
  connect(ui->checkBoxConstrain, SIGNAL(toggled(bool)), bp, SLOT(SetDrawConnectedOnly(bool)));
  connect(ui->checkBoxDrawRange, SIGNAL(toggled(bool)), bp, SLOT(SetDrawRangeEnabled(bool)));
  connect(ui->checkBoxExcludeRange, SIGNAL(toggled(bool)), bp, SLOT(SetExcludeRangeEnabled(bool)));
  connect(ui->checkBoxFill3D, SIGNAL(toggled(bool)), bp, SLOT(SetFill3D(bool)));
  connect(ui->checkBoxEraseRange, SIGNAL(toggled(bool)), bp, SLOT(SetEraseRangeEnabled(bool)));
  connect(ui->checkBoxEraseExcludeRange, SIGNAL(toggled(bool)), bp, SLOT(SetEraseExcludeRangeEnabled(bool)));
  connect(mainwnd->GetLayerCollection("MRI"), SIGNAL(LayerAdded(Layer*)), this, SLOT(UpdateWidgets()));
  connect(mainwnd->GetLayerCollection("MRI"), SIGNAL(LayerRemoved(Layer*)), this, SLOT(UpdateWidgets()));
  connect(bp, SIGNAL(FillValueChanged(double)), this, SLOT(UpdateWidgets()));
  connect(bp, SIGNAL(EraseValueChanged(double)), this, SLOT(UpdateWidgets()));
  connect(ui->pushButtonGeoClear, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegClear()));
  connect(ui->pushButtonGeoClearFilling, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegClearFilling()));
  connect(ui->pushButtonGeoGo, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegGo()));
  connect(ui->pushButtonGeoApply, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegApply()));
  connect(ui->pushButtonGeoUndo, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegUndo()));
  connect(ui->colorPickerGeoInside, SIGNAL(colorChanged(QColor)), SLOT(OnColorPickerGeoSeg(QColor)));
  connect(ui->colorPickerGeoOutside, SIGNAL(colorChanged(QColor)), SLOT(OnColorPickerGeoSeg(QColor)));
  connect(ui->colorPickerGeoFill, SIGNAL(colorChanged(QColor)), SLOT(OnColorPickerGeoSeg(QColor)));
  connect(ui->sliderGeoOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderGeoOpacity(int)));
  connect(ui->pushButtonAbort, SIGNAL(clicked(bool)), SLOT(OnButtonGeoSegAbort()));
  connect(mainwnd, SIGNAL(SupplementLayerChanged()), this, SLOT(UpdateWidgets()));

  connect(ui->pushButtonCloneCopy, SIGNAL(clicked()), mainwnd->ui->actionCopy, SLOT(trigger()));
  connect(ui->pushButtonCloneCopyStructure, SIGNAL(clicked()), mainwnd->ui->actionCopyStructure, SLOT(trigger()));
  connect(ui->pushButtonClonePaste, SIGNAL(clicked()), mainwnd->ui->actionPaste, SLOT(trigger()));

  connect(mainwnd->ui->actionCopy, SIGNAL(triggered(bool)), SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect(mainwnd->ui->actionCopyStructure, SIGNAL(triggered(bool)), SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect(mainwnd->ui->actionPaste, SIGNAL(triggered(bool)), SLOT(UpdateWidgets()), Qt::QueuedConnection);

  for (int i = 0; i < 3; i++)
  {
    RenderView2D* view = (RenderView2D*)mainwnd->GetRenderView(i);
    connect(ui->colorPickerContour, SIGNAL(colorChanged(QColor)),
            view->GetContour2D(), SLOT(SetColor(QColor)));
    connect(ui->checkBoxSmooth, SIGNAL(toggled(bool)),
            view->GetContour2D(), SLOT(SetSmooth(bool)));
    connect(view->GetContour2D(), SIGNAL(ValueChanged()), this, SLOT(UpdateWidgets()));
  }

  m_widgetsBrushSize << ui->labelBrushSize << ui->spinBoxBrushSize;

  m_widgetsReference << ui->labelReference << ui->comboBoxReference;

  m_widgetsTolerance << ui->labelTolerance << ui->spinBoxTolerance;

  //  m_widgetsConstrain << ui->checkBoxConstrain
  //                     << ui->checkBoxDrawRange
  //                     << ui->labelDrawRangeLow
  //                     << ui->labelDrawRangeHigh
  //                     << ui->lineEditDrawRangeLow
  //                     << ui->lineEditDrawRangeHigh
  //                     << ui->checkBoxExcludeRange
  //                     << ui->labelExcludeRangeHigh
  //                     << ui->labelExcludeRangeLow
  //                     << ui->lineEditExcludeRangeHigh
  //                     << ui->lineEditExcludeRangeLow;
  m_widgetsConstrain << ui->tabWidgetConstrains;

  m_widgetsSmooth << ui->checkBoxSmooth
                  << ui->labelSD
                  << ui->lineEditSD;

  m_widgetsContour << ui->labelContourColor
                   << ui->labelContourValue
                   << ui->lineEditContourValue
                   << ui->colorPickerContour
                   << ui->labelTipsContour;

  m_widgetsGeoSeg  << ui->labelGeoLambda
                   << ui->labelGeoMaxDistance
                   << ui->checkBoxMaxForegroundDistance
                   << ui->lineEditGeoMaxForegroundDistance
                   << ui->labelGeoWsize
                   << ui->lineEditGeoLambda
                   << ui->spinBoxGeoWsize
                   << ui->lineEditGeoMaxDistance
                   << ui->pushButtonGeoGo
                   << ui->pushButtonGeoClear
                   << ui->pushButtonGeoClearFilling
                   << ui->widgetGeoColors
                   << ui->sliderGeoOpacity
                   << ui->pushButtonGeoApply
                   << ui->pushButtonGeoUndo
                   << ui->labelTipsGeoS
                   << ui->widgetBusyIndicator
                   << ui->checkBoxApplySmoothing
                   << ui->lineEditSmoothingStd
                   << ui->pushButtonAbort
                   << ui->checkBoxGeoSegOverwrite
                   << ui->labelGeoMessage;

  QTimer* timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()), this, SLOT(OnIdle()) );
  timer->start( 50 );

#ifdef Q_OS_MAC
  ui->labelTips->setText(ui->labelTips->text().replace("Ctrl +", "Cmd +"));
//  ui->labelTipsContour->setText(ui->labelTipsContour->text().replace("Ctrl +", "Cmd +"));
  if (MacHelper::IsDarkMode())
  {
      ui->actionFreeHand->setIcon(MacHelper::InvertIcon(ui->actionFreeHand->icon(), QSize(), true));
      ui->actionPolyLine->setIcon(MacHelper::InvertIcon(ui->actionPolyLine->icon(), QSize(), true));
  }
#endif

  m_bToUpdateWidgets = true;
}

ToolWindowEdit::~ToolWindowEdit()
{
  QSettings settings;
  settings.setValue("ToolWindowVoxelEdit/Position", pos()-this->parentWidget()->pos());

  delete ui;
}

void ToolWindowEdit::showEvent(QShowEvent* event)
{
  Q_UNUSED(event);
  static bool bFirstTime = true;
  if ( bFirstTime )
  {
    this->move( parentWidget()->pos() + QPoint(20,100) );
    bFirstTime = false;
  }
  UpdateReconMode();
}

void ToolWindowEdit::UpdateReconMode()
{
  MainWindow* wnd = MainWindow::GetMainWindow();
  bool bReconEdit = wnd->GetRenderView(0)->GetInteractionMode() == RenderView::IM_ReconEdit;
  ui->checkBoxReconEditing->setChecked(bReconEdit);
  this->setWindowTitle(bReconEdit ? "Recon Edit" : "Voxel Edit");
}

void ToolWindowEdit::UpdateWidgets( )
{
  m_bToUpdateWidgets = true;
}


void ToolWindowEdit::OnIdle()
{
  if ( !m_bToUpdateWidgets ) // qApp->hasPendingEvents() )
  {
    return;
  }

  m_bToUpdateWidgets = false;
  ui->labelGeoMessage->clear();
  QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( true );
  }

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  bool bReconEdit = mainwnd->GetRenderView(0)->GetInteractionMode() == RenderView::IM_ReconEdit;
  int nViewId = mainwnd->GetActiveViewId();
  if (nViewId > 2 )
  {
    nViewId = 0;
  }
  RenderView2D* view = (RenderView2D*)mainwnd->GetRenderView( nViewId );
  ui->actionColorPicker->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_ColorPicker );
  ui->actionContour->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Contour );
  ui->actionFill->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Fill );
  ui->actionLiveWire->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Livewire );
  ui->actionFreeHand->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Freehand );
  ui->actionPolyLine->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Polyline );
  ui->actionClone->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_Clone );
  ui->actionAutoSeg->setChecked( view->GetAction() == Interactor2DVoxelEdit::EM_GeoSeg );

  ui->spinBoxBrushSize->setEnabled( view->GetAction() != Interactor2DVoxelEdit::EM_Fill );
  ui->spinBoxTolerance->setEnabled( view->GetAction() == Interactor2DVoxelEdit::EM_Fill );

  BrushProperty* bp = mainwnd->GetBrushProperty();
  LayerVolumeBase* layer = bp->GetReferenceLayer();

  ui->comboBoxReference->clear();
  ui->comboBoxReference->addItem( "None" );
  LayerCollection* lc = mainwnd->GetLayerCollection( "MRI" );
  int nSel = 0;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri = (LayerMRI*)lc->GetLayer( i );
    if ( layer == mri )
    {
      nSel = i+1;
    }
    ui->comboBoxReference->addItem( mri->GetName(),  QVariant::fromValue((QObject*)mri) );
  }
  ui->comboBoxReference->setCurrentIndex( nSel );

  ChangeSpinBoxValue( ui->spinBoxBrushSize, bp->GetBrushSize() );
  ChangeSpinBoxValue( ui->spinBoxTolerance, bp->GetBrushTolerance() );
  ChangeLineEditNumber(ui->lineEditBrushValue, bp->GetFillValue());
  ChangeLineEditNumber(ui->lineEditEraseValue, bp->GetEraseValue());

  double* range = bp->GetDrawRange();
  ChangeLineEditNumber( ui->lineEditDrawRangeLow, range[0]);
  ChangeLineEditNumber( ui->lineEditDrawRangeHigh, range[1] );
  range = bp->GetExcludeRange();
  ChangeLineEditNumber( ui->lineEditExcludeRangeLow, range[0] );
  ChangeLineEditNumber( ui->lineEditExcludeRangeHigh, range[1] );

  range = bp->GetEraseRange();
  ChangeLineEditNumber( ui->lineEditEraseRangeLow, range[0]);
  ChangeLineEditNumber( ui->lineEditEraseRangeHigh, range[1] );
  range = bp->GetEraseExcludeRange();
  ChangeLineEditNumber( ui->lineEditEraseExcludeRangeLow, range[0] );
  ChangeLineEditNumber( ui->lineEditEraseExcludeRangeHigh, range[1] );

  Contour2D* c2d = view->GetContour2D();
  ChangeLineEditNumber( ui->lineEditSD, c2d->GetSmoothSD() );
  ChangeLineEditNumber( ui->lineEditContourValue, c2d->GetContourValue() );

  double* rgb = c2d->GetColor();
  ui->colorPickerContour->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );

  QVariantMap geos = bp->GetGeosSettings();
  if (geos.contains("Opacity"))
    ui->sliderGeoOpacity->setValue(geos["Opacity"].toDouble()*100);
  if (geos.contains("FillColor"))
    ui->colorPickerGeoFill->setCurrentColor(geos["FillColor"].value<QColor>());
  if (geos.contains("FillColor"))
    ui->colorPickerGeoInside->setCurrentColor(geos["ForegroundColor"].value<QColor>());
  if (geos.contains("BackgroundColor"))
    ui->colorPickerGeoOutside->setCurrentColor(geos["BackgroundColor"].value<QColor>());

  int nAction = view->GetAction();
  ShowWidgets( m_widgetsBrushSize, nAction != Interactor2DVoxelEdit::EM_Contour &&
      nAction != Interactor2DVoxelEdit::EM_ColorPicker &&
      nAction != Interactor2DVoxelEdit::EM_Fill );
  //  ShowWidgets( m_widgetsReference, nAction == Interactor2DVoxelEdit::EM_Fill ||
  //                                   nAction == Interactor2DVoxelEdit::EM_Contour );
  ShowWidgets( m_widgetsTolerance, nAction == Interactor2DVoxelEdit::EM_Fill );
  ShowWidgets( m_widgetsConstrain, nAction != Interactor2DVoxelEdit::EM_ColorPicker &&
      nAction != Interactor2DVoxelEdit::EM_Contour );
  ShowWidgets( m_widgetsSmooth, nAction == Interactor2DVoxelEdit::EM_Contour );
  ShowWidgets( m_widgetsContour, nAction == Interactor2DVoxelEdit::EM_Contour );
  ShowWidgets( m_widgetsGeoSeg, nAction == Interactor2DVoxelEdit::EM_GeoSeg);
//  ui->widgetBusyIndicator->setVisible(ui->pushButtonGeoGo->isVisible() && !ui->pushButtonGeoGo->isEnabled());
  ui->widgetBusyIndicator->hide();

  ui->labelGeoLambda->hide();
  ui->labelGeoWsize->hide();
  ui->lineEditGeoLambda->hide();
  ui->spinBoxGeoWsize->hide();

  ui->widgetClone->setVisible(nAction == Interactor2DVoxelEdit::EM_Clone);

  ui->checkBoxFill3D->setVisible(nAction != Interactor2DVoxelEdit::EM_GeoSeg && nAction != Interactor2DVoxelEdit::EM_Contour);

  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( false );
  }

  ui->checkBoxSmooth->setChecked( c2d->GetSmooth() );
  ui->checkBoxConstrain->setChecked( bp->GetDrawConnectedOnly() );
  ui->checkBoxDrawRange->setChecked( bp->GetDrawRangeEnabled() );
  ui->checkBoxExcludeRange->setChecked( bp->GetExcludeRangeEnabled() );
  ui->checkBoxEraseRange->setChecked(bp->GetEraseRangeEnabled());
  ui->checkBoxEraseExcludeRange->setChecked(bp->GetEraseExcludeRangeEnabled());

  if (bReconEdit)
  {
    ui->lineEditExcludeRangeLow->setEnabled(false);
    ui->lineEditExcludeRangeHigh->setEnabled(false);
  }

//  LayerMRI* mri_draw = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
//  ui->pushButtonGeoUndo->setEnabled(mri_draw && mri_draw->HasUndo());

  LayerMRI* mri = ( LayerMRI* )mainwnd->GetActiveLayer("MRI");
  int nWnd = mainwnd->GetActiveViewId();
  ui->pushButtonCloneCopy->setEnabled( mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3 );
  ui->pushButtonCloneCopyStructure->setEnabled(mri && mri->IsVisible() && nWnd >= 0 && nWnd < 3);
  ui->pushButtonClonePaste->setEnabled( mri && mri->IsVisible() && mri->IsEditable() &&
                               nWnd >= 0 && nWnd < 3 && mri->IsValidToPaste( nWnd ) );
}

void ToolWindowEdit::OnEditMode(QAction *act)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  mainwnd->SetAction( act->data().toInt() );
  BrushProperty* bp = mainwnd->GetBrushProperty();
  if (act->data().toInt() == Interactor2DVoxelEdit::EM_GeoSeg && bp->GetReferenceLayer() == NULL)
  {
    QList<Layer*> layers = mainwnd->GetLayers("MRI");
    foreach (Layer* layer, layers)
    {
      LayerMRI* mri = (LayerMRI*)layer;
      if (mri != mainwnd->GetActiveLayer("MRI") && mri->IsVisible() &&
          mri->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT)
      {
        bp->SetReferenceLayer(mri);
        break;
      }
    }
  }
  setWindowTitle(tr("Voxel Edit - %1").arg(act->text()));

  UpdateWidgets();
}

void ToolWindowEdit::OnLineEditContourValue(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK && value > 0 )
  {
    BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
    LayerMRI* mri = (LayerMRI*)bp->GetReferenceLayer();
    for ( int i = 0; i < 3; i++ )
    {
      RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindow()->GetRenderView( i );
      Contour2D* c2d = view->GetContour2D();
      if ( c2d->GetInputImage() )
      {
        c2d->SetContourValue( value );
      }
      else if ( mri )
      {
        c2d->SetInput( mri->GetSliceImageData( view->GetViewPlane() ), value, mri->GetSlicePosition()[i], mri->GetActiveFrame() );
        c2d->SetVisible( true );
      }
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnLineEditSmoothSD(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK && value > 0 )
  {
    for ( int i = 0; i < 3; i++ )
    {
      RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindow()->GetRenderView( i );
      Contour2D* c2d = view->GetContour2D();
      c2d->SetSmoothSD( value );
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnLineEditFillValue(const QString &strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    MainWindow::GetMainWindow()->GetBrushProperty()->SetFillValue(value);
  }
}

void ToolWindowEdit::OnLineEditEraseValue(const QString &strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    MainWindow::GetMainWindow()->GetBrushProperty()->SetEraseValue(value);
  }
}

void ToolWindowEdit::OnDrawRangeChanged(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
    double* range = bp->GetDrawRange();
    if ( sender() == ui->lineEditDrawRangeLow)
    {
      bp->SetDrawRange( value, range[1] );
    }
    else if (sender() == ui->lineEditDrawRangeHigh)
    {
      bp->SetDrawRange( range[0], value );
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnExcludeRangeChanged(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
    double* range = bp->GetExcludeRange();
    if ( sender() == ui->lineEditExcludeRangeLow)
    {
      bp->SetExcludeRange( value, range[1] );
    }
    else if (sender() == ui->lineEditExcludeRangeHigh)
    {
      bp->SetExcludeRange( range[0], value );
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnEraseRangeChanged(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
    double* range = bp->GetEraseRange();
    if ( sender() == ui->lineEditEraseRangeLow)
    {
      bp->SetEraseRange( value, range[1] );
    }
    else if (sender() == ui->lineEditEraseRangeHigh)
    {
      bp->SetEraseRange( range[0], value );
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnEraseExcludeRangeChanged(const QString& strg)
{
  bool bOK;
  double value = strg.toDouble(&bOK);
  if ( bOK )
  {
    BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
    double* range = bp->GetEraseExcludeRange();
    if ( sender() == ui->lineEditEraseExcludeRangeLow)
    {
      bp->SetEraseExcludeRange( value, range[1] );
    }
    else if (sender() == ui->lineEditEraseExcludeRangeHigh)
    {
      bp->SetEraseExcludeRange( range[0], value );
    }
    UpdateWidgets();
  }
}

void ToolWindowEdit::OnComboReference(int sel)
{
  LayerVolumeBase* layer = qobject_cast<LayerVolumeBase*>(ui->comboBoxReference->itemData(sel).value<QObject*>());
  BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
  bp->SetReferenceLayer( layer );
}

void ToolWindowEdit::OnReplaceLabel()
{
  DialogReplaceLabel dlg(this);
  if (dlg.exec() == QDialog::Accepted)
  {
    LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
    int nPlane = MainWindow::GetMainWindow()->GetMainViewId();
    if (nPlane == 3 || !dlg.ReplaceSingleSlice())
      nPlane = -1;
    if (mri)
      mri->ReplaceVoxelValue(dlg.GetOriginalValue(), dlg.GetNewValue(), nPlane);
  }
}

void ToolWindowEdit::OnCheckReconEditing(bool bRecon)
{
  static int old_erase_value = 0;
  static bool exclude_enabled = false;
  static double exclude_range[2] = {0, 0};
  BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
  if (bRecon)
  {
    /*
        QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
        foreach (Layer* layer, layers)
        {
            LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
            if (mri)
            {
                mri->SetFillValue(255);
                mri->SetBlankValue(1.0);
            }
        }
        */
    old_erase_value = bp->GetEraseValue();
    double* r = bp->GetExcludeRange();
    exclude_range[0] = r[0];
    exclude_range[1] = r[1];
    exclude_enabled = bp->GetExcludeRangeEnabled();
    bp->SetFillValue(255);
    bp->SetEraseValue(1);
    bp->SetExcludeRangeEnabled(true);
    bp->SetExcludeRange(5, 250);
  }
  else
  {
    bp->SetEraseValue(old_erase_value);
    bp->SetExcludeRangeEnabled(exclude_enabled);
    bp->SetExcludeRange(exclude_range[0], exclude_range[1]);
  }
}

void ToolWindowEdit::OnButtonGeoSegClear()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
  if (mri)
  {
    mri->SaveForUndo();
    mri->ClearVoxels();
    MainWindow::GetMainWindow()->RequestRedraw();
  }
}

void ToolWindowEdit::OnButtonGeoSegUndo()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
  if (mri)
  {
    mri->Undo();
    MainWindow::GetMainWindow()->RequestRedraw();
  }
}

void ToolWindowEdit::OnButtonGeoSegClearFilling()
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
  if (mri)
  {
    mri->ClearVoxels();
    MainWindow::GetMainWindow()->RequestRedraw();
  }
}

void ToolWindowEdit::OnButtonGeoSegGo()
{
  ui->labelGeoMessage->clear();
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
  if (mri)
  {
    LayerMRI* mri_draw = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
    LayerMRI* mri_fill = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
    if (mri_draw && mri_fill)
    {
      double lambda = ui->lineEditGeoLambda->text().trimmed().toDouble();
      double max_dist = ui->lineEditGeoMaxDistance->text().trimmed().toDouble();
      double max_foreground_dist = ui->lineEditGeoMaxForegroundDistance->text().trimmed().toDouble();
      if (!ui->checkBoxMaxForegroundDistance->isChecked())
        max_foreground_dist = 0;
      int wsize = ui->spinBoxGeoWsize->value();
      mri_fill->ClearVoxels();
      bool ok = false;
      double std = 0;
      if (ui->checkBoxApplySmoothing->isChecked())
        std = ui->lineEditSmoothingStd->text().toDouble(&ok);
      mri_fill->GeodesicSegmentation(mri_draw, lambda, wsize, max_dist, ok?std:0,
                                     ui->checkBoxGeoSegOverwrite->isChecked() ? NULL : ((LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI")),
                                     max_foreground_dist);
      connect(mri_fill, SIGNAL(GeodesicSegmentationFinished(double)), this, SLOT(OnGeoSegFinished(double)), Qt::UniqueConnection);
      connect(mri_fill, SIGNAL(GeodesicSegmentationProgress(double)), this, SLOT(OnGeoSegProgress(double)), Qt::UniqueConnection);
      ui->pushButtonGeoGo->setEnabled(false);
      ui->pushButtonGeoApply->setEnabled(false);
//      ui->widgetBusyIndicator->show();
      ui->pushButtonAbort->setEnabled(true);
    }
    else
    {
      ui->labelGeoMessage->setText("<span style=\"color:red\">No inside pixels found</span>");
    }
  }
}

void ToolWindowEdit::OnButtonGeoSegAbort()
{
  LayerMRI* mri_fill = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
  if (mri_fill)
    mri_fill->GeodesicSegmentationAbort();
  ui->labelGeoMessage->clear();
}

void ToolWindowEdit::OnGeoSegFinished(double time_in_secs)
{
  ui->pushButtonGeoApply->setEnabled(true);
  ui->pushButtonGeoGo->setEnabled(true);
  ui->widgetBusyIndicator->hide();
  ui->pushButtonAbort->setEnabled(false);
  if (time_in_secs < 0)
  {
    LayerMRI* mri_fill = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
    ui->labelGeoMessage->setText(QString("<span style=\"color:red\">%1</span>")
                                 .arg(mri_fill?mri_fill->GetGeoSegErrorMessage():""));
  }
  else
    ui->labelGeoMessage->setText(QString::asprintf("Time taken: %.3fs", time_in_secs));
}

void ToolWindowEdit::OnGeoSegProgress(double val)
{
  ui->labelGeoMessage->setText(QString::asprintf("%.1f%%", val));
}

void ToolWindowEdit::OnButtonGeoSegApply()
{
  LayerMRI* mri_fill = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
  if (mri && mri_fill)
  {
    connect(mri, SIGNAL(GeodesicSegmentationApplied()), SLOT(OnButtonGeoSegClearFilling()), Qt::UniqueConnection);
    mri->GeodesicSegmentationApply(mri_fill);
  }
  ui->labelGeoMessage->clear();
}

void ToolWindowEdit::OnColorPickerGeoSeg(const QColor &color)
{
  BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
  if (sender() == ui->colorPickerGeoFill)
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
    if (mri)
    {
      QMap<int, QColor> colors;
      colors[0] = QColor(0,0,0,0);
      colors[1] = ui->colorPickerGeoFill->currentColor();
      mri->GetProperty()->SetCustomColors(colors);
      bp->SetGeosSettings("FillColor", colors[1]);
    }
  }
  else
  {
    LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
    if (mri)
    {
      QMap<int, QColor> colors;
      colors[0] = QColor(0,0,0,0);
      colors[1] = ui->colorPickerGeoInside->currentColor();
      colors[2] = ui->colorPickerGeoOutside->currentColor();
      mri->GetProperty()->SetCustomColors(colors);
      bp->SetGeosSettings("ForegroundColor", colors[1]);
      bp->SetGeosSettings("BackgroundColor", colors[2]);
    }
  }
  MainWindow::GetMainWindow()->RequestRedraw();
}

void ToolWindowEdit::OnSliderGeoOpacity(int nVal)
{
  LayerMRI* mri_fill = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_FILL"));
  LayerMRI* mri = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
  if (mri_fill)
    mri_fill->GetProperty()->SetOpacity(nVal/100.0);
  if (mri)
    mri->GetProperty()->SetOpacity(nVal/100.0);

  BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
  bp->SetGeosSettings("Opacity", nVal/100.0);
}
