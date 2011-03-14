#ifndef DIALOGWRITEMOVIEFRAMES_H
#define DIALOGWRITEMOVIEFRAMES_H

#include <QDialog>
#include <QTimer>

namespace Ui {
    class DialogWriteMovieFrames;
}

class RenderView;

class DialogWriteMovieFrames : public QDialog
{
    Q_OBJECT

public:
    explicit DialogWriteMovieFrames(QWidget *parent = 0);
    ~DialogWriteMovieFrames();

    void showEvent(QShowEvent *e);
    void closeEvent(QCloseEvent *e);

public slots:
    void UpdateUI();

signals:
    void Started();
    void Stopped();
    void Progress(int n);

protected slots:
    void OnWrite();
    void OnAbort();
    void OnOpen();

    void OnTimeOut();

private:
    Ui::DialogWriteMovieFrames *ui;
    RenderView* m_view;
    QTimer      m_timer;
    double      m_dAngleStepSize;
    int         m_nStartSlice;
    int         m_nSliceStepSize;
    QString     m_strOutputDir;
    int         m_nStepCount;
    int         m_nTotalSteps;
    bool        m_b3D;
};

#endif // DIALOGWRITEMOVIEFRAMES_H
