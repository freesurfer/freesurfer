#include "DialogSmoothSurface.h"
#include "ui_DialogSmoothSurface.h"
#include <QMessageBox>
#include "MainWindow.h"
#include "LayerSurface.h"
#include <QTimer>

DialogSmoothSurface::DialogSmoothSurface(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogSmoothSurface)
{
    ui->setupUi(this);
}

DialogSmoothSurface::~DialogSmoothSurface()
{
    delete ui;
}

void DialogSmoothSurface::OnApply()
{
  if (ValidateAll())
  {
    LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
    if ( !surf )
    {
      QMessageBox::warning(this, "Error", "No active surface found.");
      return;
    }
    int niters = ui->lineEditIterations->text().toInt();
    double lambda = ui->lineEditLambda->text().toDouble();
    double k_cutoff = ui->lineEditFrequencyCutoff->text().toDouble();
    surf->SmoothSurface(ui->comboBoxMethod->currentIndex(), niters, lambda, k_cutoff);

    QTimer::singleShot(0, MainWindow::GetMainWindow(), SIGNAL(SlicePositionChanged()));
  }
}

bool DialogSmoothSurface::ValidateAll()
{
  bool ok;
  int nVal = ui->lineEditIterations->text().toInt(&ok);
  if (!ok || nVal < 0)
  {
    QMessageBox::warning(this, "Error", "Number of iterations is not valid.");
    return false;
  }
  if (ui->comboBoxMethod->currentIndex() == 1)
    return true;

  double dVal = ui->lineEditLambda->text().toDouble(&ok);
  if (!ok || dVal < 0 || dVal > 1)
  {
    QMessageBox::warning(this, "Error", "Lambda value is not valid.");
    return false;
  }
  dVal = ui->lineEditFrequencyCutoff->text().toDouble(&ok);
  if (!ok)
  {
    QMessageBox::warning(this, "Error", "Frequency cutoff value is not valid.");
    return false;
  }
  return true;
}

void DialogSmoothSurface::OnMethod(int nMethod)
{
  ui->labelFrequencyCutoff->setVisible(nMethod == 0);
  ui->labelLambda->setVisible(nMethod == 0);
  ui->labelHelper->setVisible(nMethod == 0);
  ui->lineEditFrequencyCutoff->setVisible(nMethod == 0);
  ui->lineEditLambda->setVisible(nMethod == 0);
}
