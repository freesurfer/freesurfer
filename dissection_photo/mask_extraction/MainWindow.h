#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileInfoList>
#include <QProcess>
#include "WidgetImageView.h"
#include <QElapsedTimer>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class QLabel;

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

public slots:
  void OnButtonNext();
  void OnButtonPrev();
  void OnButtonCreateMask();
  void OnButtonClear();
  void LoadImage(int n);
  void ShowDialog();
  void OnSliderMaskOpacity(int val);

private slots:
  void OnProcessStarted();
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);

private:
  void UpdateIndex();

public:
  bool      m_bProfiling;
  QString   m_strPythonCmd;

private:
  Ui::MainWindow *ui;
  QFileInfoList  m_listInputFiles;
  QString m_strOutputFolder;
  int m_nIndex;

  QList<QColor> m_listStockColors;
  QString   m_strPyScriptPath;

  QProcess* m_proc;
  QString   m_strTempFolder;
  QElapsedTimer m_timer;
};
#endif // MAINWINDOW_H
