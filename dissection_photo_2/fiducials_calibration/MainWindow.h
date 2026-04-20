#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <QPoint>
#include <QList>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  void ShowDialog();

public slots:
  void OnButtonClear();
  void OnButtonGo();
  void OnButtonLoad();

private slots:
  void OnProcessStarted();
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);
  void OnCalibrationReady(const QList<QPoint>& pts);

public:
  QString m_strPythonCmd;

private:
  Ui::MainWindow *ui;
  QString m_strTempFolder;
  QString m_strPyScriptPath;
  QString m_strOutputFilename;
  QString m_strLastDir;

  QProcess* m_proc;
  QList<QPoint> m_listPoints;
};
#endif // MAINWINDOW_H
