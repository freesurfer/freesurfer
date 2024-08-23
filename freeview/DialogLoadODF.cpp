#include "DialogLoadODF.h"
#include "ui_DialogLoadODF.h"
#include "MainWindow.h"
#include <QMessageBox>
#include <QFileInfo>
#include <QFileDialog>
#include <QSettings>

DialogLoadODF::DialogLoadODF(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogLoadODF)
{
  ui->setupUi(this);

  QSettings s;
  ui->lineEditODF->setText(s.value("DialogLoadODF/OdfFile").toString());
  ui->lineEditMesh->setText(s.value("DialogLoadODF/MeshFile").toString());
  ui->lineEditVertexFile->setText(s.value("DialogLoadODF/VertexFile").toString());
  m_strLastDir = s.value("DialogLoadODF/LastDir").toString();
  ui->checkBoxDTK->setChecked(s.value("DialogLoadODF/DTK").toBool());
  ui->checkBoxPermuted->setChecked(s.value("DialogLoadODF/Permuted").toBool());
  ui->checkBoxHemisphere->setChecked(s.value("DialogLoadODF/Hemisphere").toBool());
}

DialogLoadODF::~DialogLoadODF()
{
  QSettings s;
  s.setValue("DialogLoadODF/OdfFile", ui->lineEditODF->text().trimmed());
  s.setValue("DialogLoadODF/MeshFile", ui->lineEditMesh->text().trimmed());
  s.setValue("DialogLoadODF/VertexFile", ui->lineEditVertexFile->text().trimmed());
  s.setValue("DialogLoadODF/LastDir", m_strLastDir);
  s.setValue("DialogLoadODF/DTK", ui->checkBoxDTK->isChecked());
  s.setValue("DialogLoadODF/Hemisphere", ui->checkBoxHemisphere->isChecked());
  s.setValue("DialogLoadODF/Permuted", ui->checkBoxPermuted->isChecked());

  delete ui;
}

void DialogLoadODF::OnButtonODF()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select ODF file",
        MainWindow::AutoSelectLastDir( m_strLastDir, "mri" ),
        "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditODF->setText(filename);
    ui->lineEditODF->setCursorPosition( ui->lineEditODF->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}

void DialogLoadODF::OnButtonVertex()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select vertex file",
        m_strLastDir,
        "Text files (*.txt);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditVertexFile->setText(filename);
    ui->lineEditVertexFile->setCursorPosition( ui->lineEditVertexFile->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}

void DialogLoadODF::OnButtonMesh()
{
  QString filename = QFileDialog::getOpenFileName(
        this,
        "Select mesh/face file",
        m_strLastDir,
        "Text files (*.txt);;All files (*)");
  if ( !filename.isEmpty() )
  {
    ui->lineEditMesh->setText(filename);
    ui->lineEditMesh->setCursorPosition( ui->lineEditMesh->text().size() );
    m_strLastDir = QFileInfo(filename).canonicalPath();
  }
}

void DialogLoadODF::OnCheckboxDTK(bool b)
{
  if (b)
    ui->checkBoxHemisphere->setChecked(true);
}

QString DialogLoadODF::GetOdfFile()
{
  QString fn = ui->lineEditODF->text().trimmed();
  if (!fn.isEmpty())
  {
    if (ui->checkBoxHemisphere->isChecked())
      fn += ":hemisphere=1";
    if (ui->checkBoxPermuted->isChecked())
      fn += ":permuted=1";
  }
  return fn;
}

QString DialogLoadODF::GetMeshFile()
{
  if (ui->checkBoxDTK->isChecked())
    return "";
  else
    return ui->lineEditMesh->text().trimmed();
}

QString DialogLoadODF::GetVertexFile()
{
  if (ui->checkBoxDTK->isChecked())
    return "";
  else
    return ui->lineEditVertexFile->text().trimmed();
}

void DialogLoadODF::OnOK()
{
  if (GetOdfFile().isEmpty())
    QMessageBox::warning(this, "Error", "Please select an ODF file");
  else if (!ui->checkBoxDTK->isChecked() && GetVertexFile().isEmpty())
    QMessageBox::warning(this, "Error", "Please select a vertex file");
  else if (!ui->checkBoxDTK->isChecked() && GetMeshFile().isEmpty())
    QMessageBox::warning(this, "Error", "Please select a mesh file");
  else
    accept();
}
