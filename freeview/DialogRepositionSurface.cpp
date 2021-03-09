/**
 * @brief Dialog window to execute surface reposition.
 *
 */
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



#include "DialogRepositionSurface.h"
#include "ui_DialogRepositionSurface.h"
#include "MainWindow.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include "LayerMRI.h"
#include <QMessageBox>
#include <QTimer>

DialogRepositionSurface::DialogRepositionSurface(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogRepositionSurface)
{
  ui->setupUi(this);
  ui->lineEditTargetX->hide();
  ui->lineEditTargetY->hide();
  ui->lineEditTargetZ->hide();
}

DialogRepositionSurface::~DialogRepositionSurface()
{
  delete ui;
}

void DialogRepositionSurface::OnApply()
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer("MRI");
  QString msg;
  if ( !surf )
    msg = "No active surface found.";
  else if ( !mri && ui->tabWidget->currentIndex() == 0)
    msg = "No active volume found.";

  if (!msg.isEmpty())
  {
    QMessageBox::warning(this, "Error", msg);
    return;
  }

  if (ValidateAll())
  {
    ui->pushButtonApply->setDisabled(true);
    if (ui->tabWidget->currentIndex() == 0)
    {
      if (ui->comboBoxTarget->currentIndex() == 0)
      {
        surf->RepositionSurface(mri, GetVertex(),
                                GetIntensity(),
                                GetNeighborSize(),
                                GetSigma(),
                                GetFlags());
      }
      else
      {
        double pos[3];
        GetCoordinate(pos);
        surf->RepositionSurface( mri, GetVertex(),
                                 pos,
                                 GetNeighborSize(),
                                 GetSigma(),
                                 GetFlags());
      }
    }
    else if (ui->tabWidget->currentIndex() == 1)
    {
      double pos[3];
      GetCoordinate(pos);
      surf->RepositionVertex(GetVertex(), pos);
    }
    else
    {
      surf->RepositionSmoothSurface(GetVertex(), GetNeighborSize(), GetSmoothingSteps());
    }
    UpdateUI();
    ui->pushButtonApply->setDisabled(false);
    QTimer::singleShot(0, MainWindow::GetMainWindow(), SIGNAL(SlicePositionChanged()));
  }
}

void DialogRepositionSurface::OnComboTarget( int nSel )
{
  //
  ui->lineEditTarget->setVisible(nSel == 0);
  ui->lineEditTargetX->setVisible(nSel == 1);
  ui->lineEditTargetY->setVisible(nSel == 1);
  ui->lineEditTargetZ->setVisible(nSel == 1);
  //  ui->labelHelper->setVisible(nSel == 0);
}

void DialogRepositionSurface::UpdateUI()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerSurface* surf = (LayerSurface*)mainwnd->GetActiveLayer( "Surface" );
  ui->pushButtonSave->setEnabled( surf && surf->IsModified() && !mainwnd->IsBusy() );
  ui->pushButtonSaveAs->setEnabled( surf && !mainwnd->IsBusy() );
  ui->pushButtonUndo->setEnabled( surf && surf->HasUndo() );
}

void DialogRepositionSurface::OnSave()
{
  MainWindow::GetMainWindow()->SaveSurface();
  UpdateUI();
}

void DialogRepositionSurface::OnSaveAs()
{
  MainWindow::GetMainWindow()->SaveSurfaceAs();
  UpdateUI();
}

int DialogRepositionSurface::GetFlags()
{
  int flags = 0;
  if (ui->radioButtonForceDirectionIn->isChecked())
    flags |= IPFLAG_FORCE_GRADIENT_IN;
  else if (ui->radioButtonForceDirectionOut->isChecked())
    flags |= IPFLAG_FORCE_GRADIENT_OUT;

  if (ui->checkBoxAllowIntersection->isChecked())
    flags |= IPFLAG_NO_SELF_INT_TEST;

  return flags;
}

int DialogRepositionSurface::GetVertex()
{
  switch (ui->tabWidget->currentIndex())
  {
  case 0:
    return ui->lineEditVertex->text().toInt();
  case 1:
    return ui->lineEditVertex2->text().toInt();
  case 2:
    return ui->lineEditVertex3->text().toInt();
  }
  return 0;
}

int DialogRepositionSurface::GetNeighborSize()
{
  if (ui->tabWidget->currentIndex() == 0)
    return ui->lineEditSize->text().toInt();
  else
    return ui->lineEditNbhdSize->text().toInt();
}

double DialogRepositionSurface::GetIntensity()
{
  return ui->lineEditTarget->text().toDouble();
}

int DialogRepositionSurface::GetSmoothingSteps()
{
  return ui->lineEditSmoothingSteps->text().toInt();
}

void DialogRepositionSurface::GetCoordinate( double* pos )
{
  if (ui->tabWidget->currentIndex() == 0)
  {
    pos[0] = ui->lineEditTargetX->text().toDouble();
    pos[1] = ui->lineEditTargetY->text().toDouble();
    pos[2] = ui->lineEditTargetZ->text().toDouble();
  }
  else
  {
    pos[0] = ui->lineEditCoordX->text().toDouble();
    pos[1] = ui->lineEditCoordY->text().toDouble();
    pos[2] = ui->lineEditCoordZ->text().toDouble();
    if (ui->radioButtonCoordRAS->isChecked())
    {
      LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
      if ( surf )
      {
        surf->GetSurfaceRASAtVertex(GetVertex(), pos);
      }
    }
  }
}

double DialogRepositionSurface::GetSigma()
{
  return ui->lineEditSigma->text().toDouble();
}

void DialogRepositionSurface::OnUndo( )
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if ( surf )
  {
    surf->Undo();
    UpdateUI();
  }
}

bool DialogRepositionSurface::ValidateAll()
{
  QString name;
  bool ok;
  if (ui->tabWidget->currentIndex() == 0)
  {
    ui->lineEditVertex->text().toInt(&ok);
    if (!ok)
      name = "Vertex";
    ui->lineEditSize->text().toInt(&ok);
    if (!ok)
      name = "Size";
    ui->lineEditSigma->text().toDouble(&ok);
    if (!ok)
      name = "Sigma";
    ui->lineEditTarget->text().toDouble(&ok);
    if ( ui->comboBoxTarget->currentIndex() == 0 && !ok )
      name = "Intensity";
    if ( ui->comboBoxTarget->currentIndex() == 1 )
    {
      ui->lineEditTargetX->text().toDouble(&ok);
      if (!ok)
        name = "Coordinate";
      ui->lineEditTargetY->text().toDouble(&ok);
      if (!ok)
        name = "Coordinate";
      ui->lineEditTargetZ->text().toDouble(&ok);
      if (!ok)
        name = "Coordinate";
    }
  }
  else if (ui->tabWidget->currentIndex() == 1)
  {
    ui->lineEditVertex2->text().toInt(&ok);
    if (!ok)
      name = "Vertex";
    ui->lineEditCoordX->text().toDouble(&ok);
    if (!ok)
      name = "Coordinate";
    ui->lineEditCoordY->text().toDouble(&ok);
    if (!ok)
      name = "Coordinate";
    ui->lineEditCoordZ->text().toDouble(&ok);
    if (!ok)
      name = "Coordinate";
  }
  else
  {
    ui->lineEditVertex3->text().toInt(&ok);
    if (!ok)
      name = "Vertex";
    ui->lineEditNbhdSize->text().toInt(&ok);
    if (!ok)
      name = "Neighborhood Size";
    ui->lineEditSmoothingSteps->text().toInt(&ok);
    if (!ok)
      name = "Smoothing Steps";
  }

  if ( !name.isEmpty() )
  {
    QMessageBox::warning(this, "Error", QString("Invalid input for ") + name);
    return false;
  }
  else
    return true;
}

void DialogRepositionSurface::OnSurfaceVertexClicked()
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if ( surf )
  {
    int nVertex = surf->GetVertexIndexAtTarget( surf->GetSlicePosition(), NULL );
    if ( nVertex >= 0 )
    {
      ui->lineEditVertex->setText( QString::number(nVertex) );
      ui->lineEditVertex2->setText(QString::number(nVertex) );
      ui->lineEditVertex3->setText(QString::number(nVertex) );
      OnCoordinateTypeChanged();
    }
  }
}

void DialogRepositionSurface::UpdateVertex()
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if (surf)
  {
    int nVertex = surf->GetVertexIndexAtTarget( surf->GetSlicePosition(), NULL );
    if ( nVertex >= 0 )
    {
      ui->lineEditVertex->setText( QString::number(nVertex) );
      ui->lineEditVertex2->setText(QString::number(nVertex) );
      ui->lineEditVertex3->setText(QString::number(nVertex) );
      OnCoordinateTypeChanged();
    }
  }
}

void DialogRepositionSurface::UpdateIntensity()
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetTopVisibleLayer("MRI");
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if (mri && surf)
  {
    double ras[3], surf_ras[3];
    mri->GetSlicePosition(ras);
    surf->GetSurfaceRASAtTarget(ras, surf_ras);
    mri->RemapPositionToRealRAS(ras, ras);
    double val = mri->GetSampledVoxelValueByRAS(ras);
    if (val >= 0)
    {
      ui->lineEditTarget->setText(QString::number(val, 'f', 2));
      OnCoordinateTypeChanged();
    }
    ui->lineEditTargetX->setText(QString::number(surf_ras[0], 'f', 2));
    ui->lineEditTargetY->setText(QString::number(surf_ras[1], 'f', 2));
    ui->lineEditTargetZ->setText(QString::number(surf_ras[2], 'f', 2));
  }
}

void DialogRepositionSurface::OnCoordinateTypeChanged()
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if ( surf )
  {
    int nVertex = ui->lineEditVertex2->text().toInt();
    if (nVertex >= 0)
    {
      double pt[3];
      if (ui->radioButtonCoordRAS->isChecked())
        surf->GetRASAtVertex(nVertex, pt);
      else
        surf->GetSurfaceRASAtVertex(nVertex, pt);
      ui->lineEditCoordX->setText(QString::number(pt[0], 'f', 2));
      ui->lineEditCoordY->setText(QString::number(pt[1], 'f', 2));
      ui->lineEditCoordZ->setText(QString::number(pt[2], 'f', 2));
    }
  }
}
