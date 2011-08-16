#ifndef WINDOWTIMECOURSE_H
#define WINDOWTIMECOURSE_H

#include <QWidget>

namespace Ui {
    class WindowTimeCourse;
}

class WindowTimeCourse : public QWidget
{
    Q_OBJECT

public:
    explicit WindowTimeCourse(QWidget *parent = 0);
    ~WindowTimeCourse();

public slots:
    void UpdateData();
    void OnFrameChanged(int n);
    void SetCurrentFrame(int n);

signals:
    void FrameChanged(int frame);

private:
    Ui::WindowTimeCourse *ui;
};

#endif // WINDOWTIMECOURSE_H
