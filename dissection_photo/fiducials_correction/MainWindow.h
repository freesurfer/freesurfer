#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  void showEvent(QShowEvent* e);

public slots:
  void UpdateWidgets();
  void OnButtonInputFolder();
  void OnButtonOutputFolder();
  void OnButtonCalibrationFile();
  void OnButtonProcess();

  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessFinished();

public:
  QString m_strPythonCmd;

private:
  void SetupScriptPath();

  Ui::MainWindow *ui;

  QString m_strPyScriptPath;
  QProcess* m_proc;
};
#endif // MAINWINDOW_H
