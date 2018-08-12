#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>

class QPushButton;
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow();

    void loadFile(const QString &fileName);

protected:
    void closeEvent(QCloseEvent *event) override;

private slots:
    void newFile();
    void open();
    bool save();
    bool saveAs();

private:
    void createActions();
};

#endif // MAINWINDOW_H
