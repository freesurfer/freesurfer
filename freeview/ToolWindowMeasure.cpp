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
#include "ToolWindowMeasure.h"
#include "ui_ToolWindowMeasure.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include <QTimer>
#include "Region2D.h"
#include "Region2DLine.h"
#include "Region2DPolyline.h"
#include "Region2DRectangle.h"
#include "SurfaceRegion.h"
#include "SurfaceRegionGroups.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "Interactor.h"
#include "Region3D.h"
#include <QSettings>
#include <QFileDialog>
#include <QMessageBox>
#include <QClipboard>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif

ToolWindowMeasure::ToolWindowMeasure(QWidget *parent) :
  QWidget(parent),
  UIUpdateHelper(),
  ui(new Ui::ToolWindowMeasure)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );

  QActionGroup* actGroup = new QActionGroup( this );
  actGroup->addAction( ui->actionContour );
  actGroup->addAction( ui->actionLabel );
  actGroup->addAction( ui->actionLine );
  actGroup->addAction( ui->actionPolyLine );
  actGroup->addAction( ui->actionRectangle );
  actGroup->addAction( ui->actionSpline );
  actGroup->addAction( ui->actionDrawOnContour );
  ui->actionContour->setVisible(false);
  ui->actionContour ->setData( Interactor::MM_SurfaceRegion );
  ui->actionLabel   ->setData( Interactor::MM_Label );
  ui->actionLine    ->setData( Interactor::MM_Line );
  ui->actionPolyLine->setData( Interactor::MM_Polyline );
  ui->actionRectangle->setData( Interactor::MM_Rectangle );
  ui->actionSpline->setData( Interactor::MM_Spline);
  ui->actionDrawOnContour->setData( Interactor::MM_DrawOnSurface );
  actGroup->setExclusive( true );
  connect(actGroup, SIGNAL(triggered(QAction*)), this, SLOT(OnAction(QAction*)) );

  m_widgets2D << ui->pushButtonCopy << ui->pushButtonExport;

  m_widgets3DDraw << ui->pushButtonSave
                 << ui->pushButtonLoad
                 << ui->colorPickerGroup
                 << ui->lineSeparator
                 << ui->pushButtonDeleteDrawing
                 << ui->pushButtonDeleteAllDrawing;

  m_widgets3D << ui->pushButtonSaveAll
              << ui->spinBoxId
              << ui->spinBoxGroup
              << ui->labelId
              << ui->labelGroup;

  foreach (QWidget* w, m_widgets3D)
    w->hide();

  m_region = NULL;
  m_surfaceRegion = NULL;
  m_3DRegion = NULL;
  m_bToUpdateWidgets = true;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  for ( int i = 0; i < 4; i++ )
  {
    RenderView* view = mainwnd->GetRenderView(i);
    if (i < 3)
    {
      connect( ((RenderView2D*)view), SIGNAL(RegionRemoved(Region2D*)),
               this, SLOT(SetRegion()));
      connect( ((RenderView2D*)view), SIGNAL(RegionSelected(Region2D*)),
               this, SLOT(SetRegion(Region2D*)));
    }
    else
    {
      connect( ((RenderView3D*)view), SIGNAL(SurfaceRegionSelected(SurfaceRegion*)),
               this, SLOT(SetSurfaceRegion(SurfaceRegion*)));
      connect( ((RenderView3D*)view), SIGNAL(SurfaceRegionRemoved(SurfaceRegion*)),
               this, SLOT(SetSurfaceRegion()));
      connect( ((RenderView3D*)view), SIGNAL(Region3DSelected(Region3D*)),
               this, SLOT(Set3DRegion(Region3D*)));
      connect( ((RenderView3D*)view), SIGNAL(Region3DRemoved(Region3D*)),
               this, SLOT(Set3DRegion()));
    }
  }
  LayerCollection* col_mri = mainwnd->GetLayerCollection("MRI");
  connect( col_mri, SIGNAL(LayerAdded(Layer*)), this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect( col_mri, SIGNAL(LayerRemoved(Layer*)), this, SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect( col_mri, SIGNAL(LayerPropertyChanged()), this, SLOT(UpdateWidgets()), Qt::QueuedConnection);

  QTimer* timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()), this, SLOT(OnIdle()) );
  timer->start( 50 );

#ifdef Q_OS_MAC
  if (MacHelper::IsDarkMode())
  {
      ui->actionLine->setIcon(MacHelper::InvertIcon(ui->actionLine->icon(), QSize(), true));
      ui->actionPolyLine->setIcon(MacHelper::InvertIcon(ui->actionPolyLine->icon(), QSize(), true));
      ui->actionSpline->setIcon(MacHelper::InvertIcon(ui->actionSpline->icon(), QSize(), true));
      ui->actionContour->setIcon(MacHelper::InvertIcon(ui->actionContour->icon(), QSize(), true));
      ui->actionDrawOnContour->setIcon(MacHelper::InvertIcon(ui->actionDrawOnContour->icon(), QSize(), true));
  }
#endif
}

ToolWindowMeasure::~ToolWindowMeasure()
{
  QSettings settings;
  settings.setValue("ToolWindowMeasure/Position", pos()-this->parentWidget()->pos());

  delete ui;
}

void ToolWindowMeasure::showEvent(QShowEvent* event)
{
  Q_UNUSED(event);
  static bool bFirstTime = true;
  if ( bFirstTime )
  {
    this->move( parentWidget()->pos() + QPoint(20,100) );
    bFirstTime = false;
    this->resize(sizeHint());
  }
}

void ToolWindowMeasure::OnAction(QAction *act)
{
  MainWindow::GetMainWindow()->SetAction( act->data().toInt() );
  setWindowTitle(tr("Measure - %1").arg(act->text()));
  UpdateWidgets();
}

void ToolWindowMeasure::SetRegion( Region2D* reg )
{
  if ( m_region )
  {
    m_region->disconnect( this );
  }
  m_region = reg;
  if ( m_region )
  {
    connect( m_region, SIGNAL(StatsUpdated()), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
    if ( m_surfaceRegion )
    {
      m_surfaceRegion->disconnect( this );
      m_surfaceRegion = NULL;
    }
    if ( m_3DRegion )
    {
      m_3DRegion->disconnect( this );
      m_3DRegion = NULL;
    }
    RenderView* view = MainWindow::GetMainWindow()->GetRenderView( 0 );
    if ( view->GetAction() == Interactor::MM_SurfaceRegion )
    {
      MainWindow::GetMainWindow()->SetAction( Interactor::MM_Line );
    }
  }
  UpdateWidgets();
}

void ToolWindowMeasure::SetSurfaceRegion( SurfaceRegion* reg )
{
  m_surfaceRegion = reg;
  if ( m_surfaceRegion )
  {
    connect(m_surfaceRegion, SIGNAL(ColorChanged(QColor)), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
    if ( m_region )
    {
      m_region->disconnect( this );
      m_region = NULL;
    }
    MainWindow::GetMainWindow()->SetAction( Interactor::MM_SurfaceRegion );
  }
  UpdateWidgets();
}

void ToolWindowMeasure::Set3DRegion( Region3D* reg )
{
  m_3DRegion = reg;
  if ( m_3DRegion )
  {
    connect(m_3DRegion, SIGNAL(ColorChanged(QColor)), this, SLOT(UpdateWidgets()), Qt::UniqueConnection);
    if ( m_region )
    {
      m_region->disconnect( this );
      m_region = NULL;
    }
    MainWindow::GetMainWindow()->SetAction( Interactor::MM_DrawOnSurface );
  }
  UpdateWidgets();
}

void ToolWindowMeasure::UpdateWidgets( )
{
  m_bToUpdateWidgets = true;
}

QString ToolWindowMeasure::GetLabelStats()
{
  QString strg;
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  LayerMRI* label = NULL, *mri = NULL;
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)lc->GetLayer( i ) )->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT )
    {
      label = ( (LayerMRI*)lc->GetLayer( i ) );
      break;
    }
  }
  for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)lc->GetLayer( i ) )->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT )
    {
      mri = ( (LayerMRI*)lc->GetLayer( i ) );
      break;
    }
  }
  if ( label && mri )
  {
    int nPlane = MainWindow::GetMainWindow()->GetMainViewId();
    if ( nPlane < 3 )
    {
      std::vector<int> ids, numbers;
      std::vector<double> means, sds;
      mri->GetLabelStats( label, nPlane, ids, numbers, means, sds );
      strg = "Id \tCount \tMean \t+/-SD\n";
      for ( size_t i = 0; i < ids.size(); i++ )
      {
        QString snum = QString("%1").arg(numbers[i], -4);
        QString smean = QString("%1").arg(means[i], -4);
        strg += QString("%1 \t%2 \t%3 \t%4\n").arg(ids[i]).arg(snum).arg(smean).arg(sds[i]);
      }
    }
  }

  return strg;
}

void ToolWindowMeasure::OnIdle()
{
  if ( !m_bToUpdateWidgets || false ) // qApp->hasPendingEvents() )
  {
    return;
  }

  // update all widgets and actions
  QList<QWidget*> allwidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( true );
  }
  RenderView* view = MainWindow::GetMainWindow()->GetRenderView( 0 );
  int nAction = view->GetAction();
  ui->actionLine->setChecked( nAction == Interactor::MM_Line );
  ui->actionPolyLine->setChecked( nAction == Interactor::MM_Polyline );
  ui->actionSpline->setChecked( nAction == Interactor::MM_Spline );
  ui->actionRectangle->setChecked( nAction == Interactor::MM_Rectangle );
  ui->actionLabel->setChecked( nAction == Interactor::MM_Label );
  ui->actionContour->setChecked( nAction == Interactor::MM_SurfaceRegion );
  ui->actionDrawOnContour->setChecked( nAction == Interactor::MM_DrawOnSurface );
  bool bLabelExist = false;
  LayerCollection* col_mri = MainWindow::GetMainWindow()->GetLayerCollection("MRI");
  for ( int i = 0; i < col_mri->GetNumberOfLayers(); i++ )
  {
    if ( ( (LayerMRI*)col_mri->GetLayer( i ) )->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT )
    {
      bLabelExist = true;
      break;
    }
  }
  ui->actionLabel->setEnabled( MainWindow::GetMainWindow()->GetMainViewId() < 3 &&
                               col_mri->GetNumberOfLayers() > 1 && bLabelExist );

  QString strg;
  RenderView2D* view2d = qobject_cast<RenderView2D*>(MainWindow::GetMainWindow()->GetMainView());
  QList<Region2D*> regions;
  if (view2d)
    regions = view2d->GetRegions();
  if ( nAction == Interactor::MM_Label )
  {
    strg = GetLabelStats();
  }
  else if (!m_region || !qobject_cast<Region2DRectangle*>(m_region))
  {
    foreach (Region2D* r, regions)
    {
      strg += r->GetShortStats() + "\n";
    }
  }
  else if ( m_region )
  {
    QStringList strgs = m_region->GetLongStats();
    for ( int i = 0; i < strgs.size(); i++ )
    {
      strg += strgs[i] + "\n";
    }
  }
  ui->textBrowserInfo->setText( strg );
  ui->pushButtonCopy->setEnabled( !strg.isEmpty() );
  ui->pushButtonExport->setEnabled( !strg.isEmpty() );
  ui->pushButtonUpdate->setVisible( nAction == Interactor::MM_Label );

  ShowWidgets( m_widgets3DDraw, nAction == Interactor::MM_DrawOnSurface || nAction == Interactor::MM_SurfaceRegion );
  ShowWidgets( m_widgets3D, nAction == Interactor::MM_SurfaceRegion);
  ShowWidgets( m_widgets2D, nAction != Interactor::MM_SurfaceRegion && nAction != Interactor::MM_DrawOnSurface);

  if ( m_surfaceRegion && ui->actionContour->isChecked() )
  {
    ui->spinBoxId->setValue( m_surfaceRegion->GetId() );
    ui->spinBoxGroup->setValue( m_surfaceRegion->GetGroup() );
    ui->spinBoxGroup->setRange( 1,
                                m_surfaceRegion->GetMRI()->GetSurfaceRegionGroups()
                                ->GetGroupIdRange( m_surfaceRegion ) );
    ui->colorPickerGroup->setCurrentColor( m_surfaceRegion->GetColor() );
  }
  if ( m_3DRegion && ui->actionDrawOnContour->isChecked() )
  {
    ui->colorPickerGroup->setCurrentColor( m_3DRegion->GetColor() );
  }

  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  ui->actionContour->setEnabled(mri &&  mri->GetProperty()->GetShowAsContour());
  bool bSurfaceRegionValid = ( mri && mri->GetProperty()->GetShowAsContour() && (mri->GetNumberOfSurfaceRegions() > 0 || mri->GetNumberOf3DRegions() > 0) );
  if ( bSurfaceRegionValid )
  {
    ui->spinBoxId->setRange( 1, mri->GetNumberOfSurfaceRegions() );
  }

  ui->pushButtonSave->setEnabled( m_surfaceRegion && bSurfaceRegionValid );
  ui->spinBoxId->setEnabled( m_surfaceRegion && bSurfaceRegionValid );
  ui->spinBoxGroup->setEnabled( m_surfaceRegion && bSurfaceRegionValid );
  ui->colorPickerGroup->setEnabled( bSurfaceRegionValid );
  ui->pushButtonSaveAll->setEnabled( bSurfaceRegionValid );
  ui->pushButtonSave->setEnabled(bSurfaceRegionValid);

  ui->pushButtonDeleteDrawing->setEnabled(m_3DRegion);
  ui->pushButtonDeleteAllDrawing->setEnabled(mri && mri->GetNumberOf3DRegions() > 0);

  m_bToUpdateWidgets = false;

  for ( int i = 0; i < allwidgets.size(); i++ )
  {
    allwidgets[i]->blockSignals( false );
  }
}

void ToolWindowMeasure::OnLoad()
{
  QString filename = QFileDialog::getOpenFileName( this,
                                                   "Load region(s) from file",
                                                   "",
                                                   "All files (*)");

  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( mri && !filename.isEmpty() )
  {
    if (ui->actionContour->isChecked() && !mri->LoadSurfaceRegions( filename ) )
    {
      QMessageBox::warning(this, "Error", QString("Can not load file ") + filename);
    }
    else if (ui->actionDrawOnContour->isChecked() && !mri->Load3DRegions( filename ) )
    {
      QMessageBox::warning(this, "Error", QString("Can not load file ") + filename);
    }
    UpdateWidgets();
  }
}

void ToolWindowMeasure::OnSave()
{
  QString filename = QFileDialog::getSaveFileName( this,
                                                   "Save region",
                                                   "",
                                                   "All files (*)");
  if (ui->actionContour->isChecked())
  {
    if ( m_surfaceRegion && !filename.isEmpty() )
    {
      if ( !m_surfaceRegion->Write( filename ) )
      {
        QMessageBox::warning(this, "Error", QString("Can not write to file ") + filename);
      }
    }
  }
  else if (ui->actionDrawOnContour->isChecked() && !filename.isEmpty())
  {
    LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
    if ( !mri->SaveAll3DRegions( filename ) )
    {
      QMessageBox::warning(this, "Error", QString("Failed to write to file ") + filename);
    }
  }
}

void ToolWindowMeasure::OnSaveAll()
{
  QString filename = QFileDialog::getSaveFileName( this,
                                                   "Save region",
                                                   "",
                                                   "All files (*)");
  if (!filename.isEmpty() )
  {
    LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
    if ( !mri->SaveAllSurfaceRegions( filename ) )
    {
      QMessageBox::warning(this, "Error", QString("Failed to write to file ") + filename);
    }
  }
}

void ToolWindowMeasure::OnUpdate()
{
  UpdateWidgets();
}

void ToolWindowMeasure::OnCopy()
{
  QApplication::clipboard()->setText(ui->textBrowserInfo->toPlainText());
}

void ToolWindowMeasure::OnExport()
{
  QString filename = QFileDialog::getSaveFileName( this,
                                                   "Export stats to file",
                                                   "",
                                                   "All files (*)");
  if (!filename.isEmpty())
  {
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      QMessageBox::warning(this, "Error", QString("Unable to write to file ") + filename);
      return;
    }

    QTextStream out(&file);
    out << ui->textBrowserInfo->toPlainText();
  }
}

void ToolWindowMeasure::OnSpinBoxId(int val)
{
  RenderView3D* view = ( RenderView3D* )MainWindow::GetMainWindow()->GetRenderView( 3 );
  view->PickSelectRegion( val );
  view->RequestRedraw();
}

void ToolWindowMeasure::OnSpinBoxGroup(int val)
{
  if ( m_surfaceRegion )
  {
    m_surfaceRegion->SetGroup( val );
    UpdateWidgets();
  }
}

void ToolWindowMeasure::OnColorGroup( const QColor& color )
{
  if ( m_surfaceRegion && ui->actionContour->isChecked())
  {
    m_surfaceRegion->GetMRI()->GetSurfaceRegionGroups()
        ->SetGroupColor( m_surfaceRegion->GetGroup(), color );
    UpdateWidgets();
  }
  else if (m_3DRegion && ui->actionDrawOnContour->isChecked())
  {
    m_3DRegion->SetColor(color);
    UpdateWidgets();
  }
}

void ToolWindowMeasure::OnButtonDelete3DRegion()
{
  RenderView3D* view = ( RenderView3D* )MainWindow::GetMainWindow()->GetRenderView( 3 );
  view->DeleteCurrent3DRegion();
}

void ToolWindowMeasure::OnButtonDeleteAll3DRegions()
{
  RenderView3D* view = ( RenderView3D* )MainWindow::GetMainWindow()->GetRenderView( 3 );
  m_3DRegion = NULL;
  view->DeleteAll3DRegions();
  UpdateWidgets();
}
