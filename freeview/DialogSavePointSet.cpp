#include "DialogSavePointSet.h"
#include "ui_DialogSavePointSet.h"
#include <QMessageBox>
#include <QFileDialog>
#include "LayerPropertyPointSet.h"

DialogSavePointSet::DialogSavePointSet(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogSavePointSet)
{
    ui->setupUi(this);
}

DialogSavePointSet::~DialogSavePointSet()
{
    delete ui;
}

void DialogSavePointSet::SetFileName(const QString &fn)
{
    ui->lineEditFileName->setText(fn);
    ui->lineEditFileName->setCursorPosition(fn.length());
}

QString DialogSavePointSet::GetFileName()
{
    return ui->lineEditFileName->text().trimmed();
}

void DialogSavePointSet::OnOK()
{
    if (GetFileName().isEmpty())
    {
        QMessageBox::information(this, "Error", "Please enter a file name to save.");
        return;
    }
    else
        accept();
}

void DialogSavePointSet::OnOpen()
{
    QString old_fn = GetFileName();
    QString fn = QFileDialog::getSaveFileName(this, "Select File To Save",
                                              (old_fn.isEmpty()?m_strLastDir:old_fn),                                              "All Files (*.*)");
    if (!fn.isEmpty())
    {
        SetFileName(fn);
    }
}

void DialogSavePointSet::SetType(int nType)
{
    ui->radioButtonControlPoint->setChecked(nType == LayerPropertyPointSet::ControlPoint);
    ui->radioButtonWayPoint->setChecked(nType == LayerPropertyPointSet::WayPoint);
}

int DialogSavePointSet::GetType()
{
    if ( ui->radioButtonControlPoint->isChecked())
        return LayerPropertyPointSet::ControlPoint;
    else
        LayerPropertyPointSet::WayPoint;
}
