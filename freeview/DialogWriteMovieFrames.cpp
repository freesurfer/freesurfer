#include "DialogWriteMovieFrames.h"
#include "ui_DialogWriteMovieFrames.h"

DialogWriteMovieFrames::DialogWriteMovieFrames(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogWriteMovieFrames)
{
    ui->setupUi(this);
}

DialogWriteMovieFrames::~DialogWriteMovieFrames()
{
    delete ui;
}
