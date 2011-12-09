/**
 * @file  ToolWindowEdit.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/09 22:09:06 $
 *    $Revision: 1.31 $
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
#include "Contour2D.h"
#include "BrushProperty.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "RenderView2D.h"
#include "DialogReplaceLabel.h"
#include <QTimer>
#include <QSettings>
#include <QDebug>

ToolWindowEdit::ToolWindowEdit(QWidget *parent) :
  QWidget(parent),
  UIUpdateHelper(),
  ui(new Ui::ToolWindowEdit)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );
  QActionGroup* ag = new QActionGroup( this );
  ag->addAction( ui->actionContour );
  ag->addAction( ui->actionFill );
  ag->addAction( ui->actionFreeHand );
  ag->addAction( ui->actionLiveWire );
  ag->addAction( ui->actionPolyLine );
  ag->addAction( ui->actionColorPicker );
  ag->setExclusive( true );
  ui->actionContour->setData( Interactor2DVoxelEdit::EM_Contour );
  ui->actionColorPicker->setData( Interactor2DVoxelEdit::EM_ColorPicker );
  ui->actionFill->setData( Interactor2DVoxelEdit::EM_Fill );
  ui->actionFreeHand->setData( Interactor2DVoxelEdit::EM_Freehand );
  ui->actionLiveWire->setData( Interactor2DVoxelEdit::EM_Livewire );
  ui->actionPolyLine->setData( Interactor2DVoxelEdit::EM_Polyline );
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
  connect(mainwnd->GetLayerCollection("MRI"), SIGNAL(LayerAdded(Layer*)), this, SLOT(UpdateWidgets()));
  connect(mainwnd->GetLayerCollection("MRI"), SIGNAL(LayerRemoved(Layer*)), this, SLOT(UpdateWidgets()));
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

  m_widgetsConstrain << ui->checkBoxConstrain
                     << ui->checkBoxDrawRange
                     << ui->labelDrawRangeLow
                     << ui->labelDrawRangeHigh
                     << ui->lineEditDrawRangeLow
                     << ui->lineEditDrawRangeHigh
                     << ui->checkBoxExcludeRange
                     << ui->labelExcludeRangeHigh
                     << ui->labelExcludeRangeLow
                     << ui->lineEditExcludeRangeHigh
                     << ui->lineEditExcludeRangeLow;

  m_widgetsSmooth << ui->checkBoxSmooth
                  << ui->labelSD
                  << ui->lineEditSD;

  m_widgetsContour << ui->labelContourColor
                   << ui->labelContourValue
                   << ui->lineEditContourValue
                   << ui->colorPickerContour;

  QTimer* timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()), this, SLOT(OnIdle()) );
  timer->start( 50 );

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
  static bool bFirstTime = true;
  if ( bFirstTime )
  {
    QSettings settings;
    QVariant v = settings.value( "ToolWindowVoxelEdit/Position", QPoint( 200, 20 ) );
    this->move( parentWidget()->pos() + v.toPoint() );
    bFirstTime = false;
  }
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

  QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( true );
  }

  MainWindow* mainwnd = MainWindow::GetMainWindow();
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

  double* range = bp->GetDrawRange();
  ChangeLineEditNumber( ui->lineEditDrawRangeLow, range[0]);
  ChangeLineEditNumber( ui->lineEditDrawRangeHigh, range[1] );
  range = bp->GetExcludeRange();
  ChangeLineEditNumber( ui->lineEditExcludeRangeLow, range[0] );
  ChangeLineEditNumber( ui->lineEditExcludeRangeHigh, range[1] );

  Contour2D* c2d = view->GetContour2D();
  ChangeLineEditNumber( ui->lineEditSD, c2d->GetSmoothSD() );
  ChangeLineEditNumber( ui->lineEditContourValue, c2d->GetContourValue() );

  double* rgb = c2d->GetColor();
  ui->colorPickerContour->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );

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

  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( false );
  }

  ui->checkBoxSmooth->setChecked( c2d->GetSmooth() );
  ui->checkBoxConstrain->setChecked( bp->GetDrawConnectedOnly() );
  ui->checkBoxDrawRange->setChecked( bp->GetDrawRangeEnabled() );
  ui->checkBoxExcludeRange->setChecked( bp->GetExcludeRangeEnabled() );

  m_bToUpdateWidgets = false;
}

void ToolWindowEdit::OnEditMode(QAction *act)
{
  MainWindow::GetMainWindow()->SetAction( act->data().toInt() );
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
