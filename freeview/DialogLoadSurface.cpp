#include "DialogLoadSurface.h"
#include "ui_DialogLoadSurface.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
#include "MainWindow.h"
#include <QSettings>

DialogLoadSurface::DialogLoadSurface(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadSurface)
{
  ui->setupUi(this);

  QSettings s;
  ui->lineEditFilename->setText(s.value("DialogLoadSurface/Filename").toString());
  ui->lineEditFilename->setCursorPosition(ui->lineEditFilename->text().size());
  ui->checkBoxInflated->setChecked(s.value("DialogLoadSurface/LoadInflated").toBool());
  ui->checkBoxOrig->setChecked(s.value("DialogLoadSurface/LoadOrig").toBool());
  ui->checkBoxPial->setChecked(s.value("DialogLoadSurface/LoadPial").toBool());
  ui->checkBoxWhite->setChecked(s.value("DialogLoadSurface/LoadWhite").toBool());
  UpdateStatus();
}

DialogLoadSurface::~DialogLoadSurface()
{
  QSettings s;
  s.setValue("DialogLoadSurface/Filename", ui->lineEditFilename->text().trimmed());
  s.setValue("DialogLoadSurface/LoadInflated", ui->checkBoxInflated->isChecked());
  s.setValue("DialogLoadSurface/LoadOrig", ui->checkBoxOrig->isChecked());
  s.setValue("DialogLoadSurface/LoadPial", ui->checkBoxPial->isChecked());
  s.setValue("DialogLoadSurface/LoadWhite", ui->checkBoxWhite->isChecked());

  delete ui;
}

void DialogLoadSurface::accept()
{
  if (ui->lineEditFilename->text().trimmed().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please enter the filename of the surface file to load");
  }
  else
    QDialog::accept();
}

void DialogLoadSurface::OnOpen()
{
  QString lastdir = QFileInfo(ui->lineEditFilename->text().trimmed()).absoluteFilePath();
  QString filename = QFileDialog::getOpenFileName( this, "Select surface file",
                                                   MainWindow::AutoSelectLastDir(lastdir, "surf" ),
                                                   "Surface files (*)");
  if (!filename.isEmpty())
  {
    ui->lineEditFilename->setText(filename);
    ui->lineEditFilename->setCursorPosition(filename.size());
  }

  UpdateStatus();
}

void DialogLoadSurface::UpdateStatus()
{
  /*
  QString fn = GetFilename();
  if (fn.contains(".inflated"))
  {
    ui->checkBoxInflated->setChecked(true);
    ui->checkBoxInflated->setEnabled(false);
  }
  else
    ui->checkBoxInflated->setEnabled(true);

  if (fn.contains(".pial"))
  {
    ui->checkBoxPial->setChecked(true);
    ui->checkBoxPial->setEnabled(false);
  }
  else
    ui->checkBoxPial->setEnabled(true);

  if (fn.contains(".white"))
  {
    ui->checkBoxWhite->setChecked(true);
    ui->checkBoxWhite->setEnabled(false);
  }
  else
    ui->checkBoxWhite->setEnabled(true);

  if (fn.contains(".orig"))
  {
    ui->checkBoxOrig->setChecked(true);
    ui->checkBoxOrig->setEnabled(false);
  }
  else
    ui->checkBoxOrig->setEnabled(true);
  */
}

QString DialogLoadSurface::GetFilename()
{
  return ui->lineEditFilename->text().trimmed();
}

QStringList DialogLoadSurface::GetSupFiles()
{
  QStringList list;
  if (ui->checkBoxInflated->isChecked())
    list << "inflated";
  if (ui->checkBoxOrig->isChecked())
    list << "orig";
  if (ui->checkBoxPial->isChecked())
    list << "pial";
  if (ui->checkBoxWhite->isChecked())
    list << "white";

  return list;
}
