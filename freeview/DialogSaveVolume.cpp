#include "DialogSaveVolume.h"
#include "ui_DialogSaveVolume.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QFile>

DialogSaveVolume::DialogSaveVolume(QWidget *parent, const QString& filepath) :
    QDialog(parent),
    ui(new Ui::DialogSaveVolume)
{
    ui->setupUi(this);
    ui->lineEditFileName->setText(filepath);
}

DialogSaveVolume::~DialogSaveVolume()
{
    delete ui;
}

QString DialogSaveVolume::GetFileName()
{
    return ui->lineEditFileName->text().trimmed();
}

bool DialogSaveVolume::GetResample()
{
    return !ui->checkBoxNoResample->isChecked();
}

void DialogSaveVolume::OnOK()
{
    if (GetFileName().isEmpty())
        QMessageBox::warning(this, "Error", "Please enter the file name to save.");
    else
    {
        if (QFile::exists(GetFileName()))
        {
            if (QMessageBox::question(this, "File exists", "File exists. Do you want to overwrite it?", QMessageBox::Yes, QMessageBox::No) ==
                QMessageBox::Yes)
                accept();
        }
        else
            accept();
    }
}

void DialogSaveVolume::OnOpen()
{
    QString fn = QFileDialog::getSaveFileName( this, "Save volume",
                      GetFileName(),
                      "Volume files (*.mgz *.mgh *.nii *.nii.gz *.img *.mnc);;All files (*.*)");
    if ( !fn.isEmpty() )
    {
      ui->lineEditFileName->setText(fn);
      ui->lineEditFileName->setCursorPosition(fn.length());
    }
}
