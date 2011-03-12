#include "DialogAbout.h"
#include "ui_DialogAbout.h"

DialogAbout::DialogAbout(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogAbout)
{
    ui->setupUi(this);
    QString strg = ui->labelVersion->text();
    strg.replace("xxx", QString(__DATE__) + " " + __TIME__);
    ui->labelVersion->setText(strg);
}

DialogAbout::~DialogAbout()
{
    delete ui;
}
