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
  void OnButtonNext();
  void OnButtonPrev();
  void OnButtonRegister();
  void OnButtonClear();
  void LoadImage(int n);
  void ShowDialog();

private slots:
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessStarted();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);

public:
  QString m_strPythonCmd;

private:
  void UpdateIndex();
  void SetupScriptPath();

  Ui::MainWindow *ui;
  QFileInfoList  m_listInputFiles;
  QString m_strOutputFolder;
  int m_nNumberOfExpectedPoints;
  int m_nIndex;
  QList< QList<QPoint> > m_listData;
  QString m_strPyScriptPath;

  QProcess* m_proc;
};
#endif // MAINWINDOW_H
