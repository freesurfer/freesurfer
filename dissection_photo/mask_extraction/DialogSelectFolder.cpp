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
}

DialogSelectFolder::~DialogSelectFolder()
{
  QSettings s;
  s.setValue("CurrentFolder/Input", GetInputPath());
  s.setValue("CurrentFolder/Output", GetOutputPath());

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

void DialogSelectFolder::OnButtonRegister()
{
  if (GetInputPath().isEmpty() || GetOutputPath().isEmpty())
  {
    QMessageBox::warning(this, "", "Please select input/output folder");
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
