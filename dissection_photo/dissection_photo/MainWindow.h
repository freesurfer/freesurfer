#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileInfoList>
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

public slots:
  void OnButtonInputFolder();
  void OnButtonOutputFolder();
  void OnButtonContinue();
  void OnButtonCalibrationFile();
  void OnTogglePointMode(bool b);
  void OnButtonPrevious();
  void OnButtonNext();
  void OnButtonProcess();
  void OnButtonClear();

private slots:
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessStarted();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);

private:
  void SetupScriptPath();
  void UpdateIndex();
  void LoadImage(int n);

  Ui::MainWindow *ui;
  QString  m_strInputFolder;
  QString  m_strOutputFolder;
  QString  m_strCalibrationFile;

  QString  m_strPythonCmd;

  QFileInfoList  m_listInputFiles;
  int m_nNumberOfExpectedPoints;
  int m_nIndex;
  QList< QList<QPoint> > m_listData;
  QString m_strPyScriptRetrospective;

  QProcess* m_proc;
};
#endif // MAINWINDOW_H
