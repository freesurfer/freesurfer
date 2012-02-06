#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
  class MainWindow;
}

class QProcess;

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

  void closeEvent(QCloseEvent *);

  enum MessageType { Normal = 0, Highlight, Error, Warning };

public slots:
  void OnButtonRun();
  bool OnButtonAbort();
  void OnButtonOpenInput();
  void OnButtonOpenOutput();
  void UpdateStatus();
  void OnStandardOutput();

protected:
  void AddMessage(const QString& text, int type = Normal);

private:
  Ui::MainWindow *ui;
  QProcess*   m_process;  // process/thread that runs mri_convert
  QString     m_strLastDir;
};

#endif // MAINWINDOW_H
