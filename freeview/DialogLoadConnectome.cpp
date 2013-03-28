#include "DialogLoadConnectome.h"
#include "ui_DialogLoadConnectome.h"

#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>

DialogLoadConnectome::DialogLoadConnectome(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogLoadConnectome)
{
    ui->setupUi(this);
    QSettings s;
    ui->lineEditCMAT->setText(s.value("DialogConnectome/CMAT").toString());
    ui->lineEditParcel->setText(s.value("DialogConnectome/Parcel").toString());
    ui->lineEditCMAT->setCursorPosition(ui->lineEditCMAT->text().size());
    ui->lineEditParcel->setCursorPosition(ui->lineEditParcel->text().size());
    QString strg = GetCMATFilename();
    if (!strg.isEmpty())
      m_strLastDir = QFileInfo(strg).absolutePath();
}

DialogLoadConnectome::~DialogLoadConnectome()
{
  QSettings s;
  s.setValue("DialogConnectome/CMAT", GetCMATFilename());
  s.setValue("DialogConnectome/Parcel", GetParcelFilename());
  s.setValue("DialogConnectome/LastDir", m_strLastDir);

  delete ui;
}

QString DialogLoadConnectome::GetCMATFilename()
{
  return ui->lineEditCMAT->text().trimmed();
}

QString DialogLoadConnectome::GetParcelFilename()
{
  return ui->lineEditParcel->text().trimmed();
}

QString DialogLoadConnectome::GetCTABFilename()
{
  return ui->lineEditCTAB->text().trimmed();
}

void DialogLoadConnectome::OnButtonOpenCMAT()
{
  QString fn = QFileDialog::getOpenFileName(this, "Select CMAT File", m_strLastDir,
                                            "CMAT files (*.cmat);;All files (*)");
  if (!fn.isEmpty())
  {
    ui->lineEditCMAT->setText(fn);
    ui->lineEditCMAT->setCursorPosition(fn.size());
    m_strLastDir = QFileInfo(fn).absolutePath();
  }
}

void DialogLoadConnectome::OnButtonOpenParcel()
{
  QString fn = QFileDialog::getOpenFileName(this, "Select Parcellation File", m_strLastDir,
                                            "Parcellation files (*.mgz);;All files (*)");
  if (!fn.isEmpty())
  {
    ui->lineEditParcel->setText(fn);
    ui->lineEditParcel->setCursorPosition(fn.size());
    m_strLastDir = QFileInfo(fn).absolutePath();
  }
}

void DialogLoadConnectome::OnButtonOpenCTAB()
{
  QString fn = QFileDialog::getOpenFileName(this, "Select Color Table File", m_strLastDir,
                                            "CTAB files (*.ctab);;All files (*)");
  if (!fn.isEmpty())
  {
    ui->lineEditCTAB->setText(fn);
    ui->lineEditCTAB->setCursorPosition(fn.size());
  }
}

void DialogLoadConnectome::OnOK()
{
  if (GetCMATFilename().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please select a CMAT file.");
    return;
  }
  if (GetParcelFilename().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please select a parcellation file.");
    return;
  }
  if (GetCTABFilename().isEmpty())
  {
    QMessageBox::warning(this, "Error", "Please select a color table file.");
    return;
  }

  accept();
}
