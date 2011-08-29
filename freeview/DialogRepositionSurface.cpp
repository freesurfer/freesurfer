/**
 * @file  DialogRepositionSurface.cpp
 * @brief Dialog window to execute surface reposition.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/08/29 15:24:59 $
 *    $Revision: 1.2 $
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



#include "DialogRepositionSurface.h"
#include "ui_DialogRepositionSurface.h"
#include "MainWindow.h"
#include "LayerSurface.h"
#include "LayerMRI.h"
#include <QMessageBox>

DialogRepositionSurface::DialogRepositionSurface(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogRepositionSurface)
{
    ui->setupUi(this);
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
  else if ( !mri )
    msg = "No active volume found.";

  if (!msg.isEmpty())
  {
    QMessageBox::warning(this, "Error", msg);
    return;
  }

  if (ValidateAll())
  {
    if (ui->comboBoxTarget->currentIndex() == 0)
    {
      surf->RepositionSurface(mri, GetVertex(),
                              GetIntensity(),
                              GetNeighborSize(),
                              GetSigma());
    }
    else
    {
      double pos[3];
      GetCoordinate(pos);
      surf->RepositionSurface( mri, GetVertex(),
                                     pos,
                                     GetNeighborSize(),
                                     GetSigma() );
    }
    UpdateUI();
  }
}

void DialogRepositionSurface::OnComboTarget( int nSel )
{
 //
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

int DialogRepositionSurface::GetVertex()
{
  return ui->lineEditVertex->text().toInt();
}

int DialogRepositionSurface::GetNeighborSize()
{
  return ui->lineEditSize->text().toInt();
}

double DialogRepositionSurface::GetIntensity()
{
  return ui->lineEditTarget->text().toDouble();
}

void DialogRepositionSurface::GetCoordinate( double* pos )
{
  QStringList list = ui->lineEditTarget->text().split(",", QString::SkipEmptyParts);
  if (list.size() < 3)
    list = ui->lineEditTarget->text().split(" ", QString::SkipEmptyParts);

  for ( int i = 0; i < 3; i++ )
    pos[i] = list[i].toDouble();
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
  long nval;
  double dval;
  bool ok;
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
    QStringList list = ui->lineEditTarget->text().split(",", QString::SkipEmptyParts);
    if (list.size() < 3)
      list = ui->lineEditTarget->text().split(" ", QString::SkipEmptyParts);
    if ( list.size() < 3 )
      name = "Coordinate";
    for (int i = 0; i < 3; i++)
    {
      list[i].toDouble(&ok);
      if (!ok)
        name = "Coordinate";
    }
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
      ui->lineEditVertex->setText( QString("%1").arg(nVertex) );
  }
}

