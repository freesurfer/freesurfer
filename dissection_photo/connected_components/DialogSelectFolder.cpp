#include "DialogSelectFolder.h"
#include "ui_DialogSelectFolder.h"
#include <QFileDialog>
#include <QSettings>
#include <QDir>
#include <QMessageBox>

DialogSelectFolder::DialogSelectFolder(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogSelectFolder)
{
  ui->setupUi(this);

  QSettings s;
  QString path = s.value("CurrentFolder/Input").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathInput->setText(path);
    ui->lineEditPathInput->setCursorPosition( path.size() );
  }
  path = s.value("CurrentFolder/Output").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathOutput->setText(path);
    ui->lineEditPathOutput->setCursorPosition( path.size() );
  }
  path = s.value("CurrentFolder/Mask").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathMask->setText(path);
    ui->lineEditPathMask->setCursorPosition( path.size() );
  }
}

DialogSelectFolder::~DialogSelectFolder()
{
  QSettings s;
  s.setValue("CurrentFolder/Input", GetInputPath());
  s.setValue("CurrentFolder/Output", GetOutputPath());
  s.setValue("CurrentFolder/Mask", GetMaskPath());

  delete ui;
}

void DialogSelectFolder::OnButtonInputPath()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", GetInputPath());
  if (!path.isEmpty())
  {
    ui->lineEditPathInput->setText(path);
    ui->lineEditPathInput->setCursorPosition( path.size() );
  }
}

void DialogSelectFolder::OnButtonOutputPath()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", GetInputPath());
  if (!path.isEmpty())
  {
    ui->lineEditPathOutput->setText(path);
    ui->lineEditPathOutput->setCursorPosition( path.size() );
  }
}

void DialogSelectFolder::OnButtonMaskPath()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", GetInputPath());
  if (!path.isEmpty())
  {
    ui->lineEditPathMask->setText(path);
    ui->lineEditPathMask->setCursorPosition( path.size() );
  }
}

void DialogSelectFolder::OnButtonRegister()
{
  if (GetInputPath().isEmpty() || GetOutputPath().isEmpty() || GetMaskPath().isEmpty())
  {
    QMessageBox::warning(this, "", "Please select input/output/mask folder");
    reject();
  }
  else if (!QDir(GetOutputPath()).exists())
  {
    QMessageBox::warning(this, "", "Please make sure output folder exist");
    reject();
  }
  else
    accept();
}

QString DialogSelectFolder::GetInputPath()
{
  return ui->lineEditPathInput->text().trimmed();
}

QString DialogSelectFolder::GetOutputPath()
{
  return ui->lineEditPathOutput->text().trimmed();
}

QString DialogSelectFolder::GetMaskPath()
{
  return ui->lineEditPathMask->text().trimmed();
}
