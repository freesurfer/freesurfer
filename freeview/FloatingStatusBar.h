#ifndef FLOATINGSTATUSBAR_H
#define FLOATINGSTATUSBAR_H

#include <QWidget>
#include <QTimer>

namespace Ui {
    class FloatingStatusBar;
}

class FloatingStatusBar : public QWidget
{
    Q_OBJECT

public:
    explicit FloatingStatusBar(QWidget *parent = 0);
    ~FloatingStatusBar();

public slots:
    void SetProgress( int nProgress );
    void ShowProgress();
    void HideProgress();
    void Reposition();

private slots:
    void OnProgressTimer();

private:
    Ui::FloatingStatusBar *ui;
    QTimer*       m_timer;

};

#endif // FLOATINGSTATUSBAR_H
