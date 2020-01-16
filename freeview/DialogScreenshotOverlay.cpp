#include "DialogScreenshotOverlay.h"
#include "ui_DialogScreenshotOverlay.h"
#include <QFileDialog>
#include <QDir>
#include "LayerSurface.h"
#include "SurfaceOverlay.h"
#include "MainWindow.h"
#include <QTimer>
#include <QMessageBox>
#include "RenderView.h"

DialogScreenshotOverlay::DialogScreenshotOverlay(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogScreenshotOverlay)
{
  ui->setupUi(this);
  connect(ui->toolButtonOpen, SIGNAL(clicked(bool)), SLOT(OnButtonOpen()));
  connect(ui->pushButtonSave, SIGNAL(clicked(bool)), SLOT(OnButtonSave()));
}

DialogScreenshotOverlay::~DialogScreenshotOverlay()
{
  delete ui;
}

void DialogScreenshotOverlay::showEvent(QShowEvent *e)
{
  ui->labelProgress->clear();
  ui->widgetControls->setEnabled(true);
  if (ui->lineEditDirectory->text().trimmed().isEmpty())
    ui->lineEditDirectory->setText(QDir::currentPath());
}

void DialogScreenshotOverlay::OnButtonOpen()
{
  QString path = QFileDialog::getExistingDirectory(this, "", ui->lineEditDirectory->text().trimmed());
  if (!path.isEmpty())
  {
    ui->lineEditDirectory->setText(path);
  }
}

void DialogScreenshotOverlay::OnButtonSave()
{
  QString format = ui->lineEditFilename->text().trimmed();
  format.replace("%s", "%1");
  if (format.isEmpty())
  {
    QMessageBox::information(this, "Error", "Please enter a valid filename format");
    return;
  }

  m_nIndex = 0;
  ui->widgetControls->setEnabled(false);
  OnTimeOut();
}

void DialogScreenshotOverlay::OnTimeOut()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  RenderView* view = MainWindow::GetMainWindow()->GetMainView();

  LayerSurface* surface = qobject_cast<LayerSurface*>(mainwnd->GetActiveLayer("Surface"));
  if (!surface || m_nIndex >= surface->GetNumberOfOverlays())
  {
    ui->widgetControls->setEnabled(true);
    return;
  }

  SettingsScreenshot settings = MainWindow::GetMainWindow()->GetScreenShotSettings();
  QString format = ui->lineEditFilename->text().trimmed();
  format.replace("%s", "%1");
  if (!format.contains("%1"))
  {
    QFileInfo fi(format);
    format = fi.completeBaseName() + ".%1." + fi.suffix();
  }
  QString fn;
  QString name = surface->GetOverlay(m_nIndex)->GetName();
  name.remove(".mgz");
  name.remove(".mgh");
  name.remove(".nii.gz");
  name.remove(".nii");
  QFileInfo fi(QDir(ui->lineEditDirectory->text()), format.arg(name));
  fn = fi.absoluteFilePath();
  surface->SetActiveOverlay(m_nIndex);
  view->RequestRedraw(true);
  view->SaveScreenShot(fn, settings.AntiAliasing, settings.Magnification, settings.AutoTrim);
  m_nIndex++;
  ui->labelProgress->setText(QString("%1/%2").arg(m_nIndex).arg(surface->GetNumberOfOverlays()));
  QTimer::singleShot(250, this, SLOT(OnTimeOut()));
}
