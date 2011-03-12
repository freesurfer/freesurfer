#ifndef DIALOGWRITEMOVIEFRAMES_H
#define DIALOGWRITEMOVIEFRAMES_H

#include <QDialog>

namespace Ui {
    class DialogWriteMovieFrames;
}

class DialogWriteMovieFrames : public QDialog
{
    Q_OBJECT

public:
    explicit DialogWriteMovieFrames(QWidget *parent = 0);
    ~DialogWriteMovieFrames();

private:
    Ui::DialogWriteMovieFrames *ui;
};

#endif // DIALOGWRITEMOVIEFRAMES_H
