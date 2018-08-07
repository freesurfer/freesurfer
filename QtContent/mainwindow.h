#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
#include <QPlainTextEdit>
#include <QFileInfo>
#include "configurationfileform.h"

QString bvecSeeker(QFileInfoList list, QFileInfo current);
QString bvalSeeker(QFileInfoList list, QFileInfo current);

class QPushButton;
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow();
    ~MainWindow();

    void loadFile(const QString &fileName);

protected:
    //void closeEvent(QCloseEvent *event) override;

private slots:
    void newFile();
    void open();
    /*bool save();
    bool saveAs();*/

private:
    void createActions();
    void createStatusBar();
    void setCurrentFile(const QString &fileName);

    QPlainTextEdit *textEdit;
    QString curFile;
    ConfigurationFileForm* table;
};

#endif // MAINWINDOW_H
