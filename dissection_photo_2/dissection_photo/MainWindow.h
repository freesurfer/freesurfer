#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileInfoList>
#include <QProcess>
#include <QVariantMap>
#include <QElapsedTimer>
#include "WidgetImageView.h"
#include "MaskProcessor.h"
#include <QFileSystemWatcher>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class QLabel;
class ProgressWindow;

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  void resizeEvent(QResizeEvent* e);
  void moveEvent(QMoveEvent* e);

  bool eventFilter(QObject *watched, QEvent *event);

public slots:
  void OnButtonInputFolder();
  void OnButtonOutputFolderCC();
  void OnButtonContinue();
  void OnButtonCalibrationFile();
  void OnTogglePointMode(bool b);
  void OnButtonPrevious();
  void OnButtonNext();
  void OnButtonProcess();
  void OnButtonClear();
  void OnButtonLoadMask();
  void OnButtonCreateMask();
  void OnLinkActivated(const QString& link);
  void OnButtonTogglePreview(bool b);
  void OnButtonCalibrationPhoto();
  void OnButtonBackToStart();
  void OnButtonSaveCalibration();

private slots:
  void OnProcessOutputMessage();
  void OnProcessErrorMessage();
  void OnProcessStarted();
  void OnProcessFinished();
  void OnProcessError(QProcess::ProcessError);
  void OnButtonProceedToSeg();
  void OnButtonProceedToCC();
  void OnButtonProceedToRecon();
  void OnButtonProceedToCalibration();
  void OnSliderSegOpacity(int n);
  void OnSliderLabelOpacity(int n);
  void OnLastRegionEdited(int n);
  void OnToggleMask();
  void OnFileChanged(const QString& path);
  void ShowProgressWindow(bool bShow = true);
  void OnThumbnailIndexChanged(int n);
  void OnButtonCalibration();
  void OnPointEdited();

private:
  void SetupScriptPath();
  void UpdateIndex();
  void LoadImage(int n, bool no_processing = false, bool bPreview = false);
  void ClearFolder(const QString& path);
  QImage NpyToImage(const QString& fn, const QSize& sz);
  QList<QPoint> GetCalibrationPointsList(const QVariantMap& info);
  void RepositionProgressWindow();
  void SetInputFolder(const QString& path, const QString& root_output = "");
  void SetCalibrationFile(const QString& path);
  void SetCalibrationPhoto(const QString& path);
  void RunRetrospectiveCorrection(const QString& infile, const QList<QPoint>& pts, double w, double h);

  Ui::MainWindow *ui;
  QString  m_strInputFolder;
  QString  m_strRootOutputFolder;
  QString  m_strOutputFolder;
  QString  m_strCalibrationFile;
  QString  m_strCalibrationPhoto;
  QString  m_strFinalOutputFolder;
  QString  m_strMaskFolder;

  QString  m_strNNUnetScriptFolder;
  QString  m_strNNUnetModelFolder;

  QString  m_strPythonCmd;

  QFileInfoList  m_listInputFiles;
  QFileInfoList  m_listOutputFiles;
  QFileInfoList  m_listMaskFiles;
  QMap<int, QList<RECT_REGION> > m_listRegionData;
  QList<QColor> m_listStockColors;

  int m_nNumberOfExpectedPoints;
  int m_nIndex;
  QMap<int, QList<QPoint> > m_listPointData;
  bool m_bCalibratedMode;
  QString m_strPyScriptRetrospective;
  QString m_strPyScriptFiducialsCorrection;
  QString m_strPyScriptFiducialsDetection;
  QString m_strPyScriptFiducialsCalibration;
  QString m_strPyScriptMaskToCC;

  QProcess* m_proc;

  QString m_sTempDir;

  MaskProcessor  m_maskProcessor;
  QElapsedTimer  m_elapsedTimer;
  QFileSystemWatcher  m_fileWatcher;
  QStringList    m_listQueuedFiles;
  QVariantMap    m_mapInputImageSize;
  ProgressWindow* m_wndProgress;
};
#endif // MAINWINDOW_H
