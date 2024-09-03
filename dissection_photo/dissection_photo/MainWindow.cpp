#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "CommonDef.h"
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <QScreen>
#include <QTemporaryDir>
#include <QTimer>
#include <QTextStream>
#include <QDebug>
#include <QDateTime>
#include "cnpy.h"
#include "ProgressWindow.h"

#define SCRIPT_RETROSPECTIVE "func_retrospective_correction.py"
#define SCRIPT_FIDUCIALS_CORRECTION "func_fiducials_correction.py"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);
  addAction(ui->actionToggleMask);
  ui->stackedWidget->setCurrentWidget(ui->pageStartUp);

  m_strPythonCmd = PY_COMMAND;

  SetupScriptPath();

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QGuiApplication::primaryScreen()->geometry().contains(rc))
    setGeometry(rc);
  ui->splitter->restoreState(s.value("MainWindow/SplitterState").toByteArray());

  QString path = s.value("CurrentFolder/Input").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathInput->setText(path);
    ui->lineEditPathInput->setCursorPosition( path.size() );
    m_strInputFolder = path;
  }
  path = s.value("CurrentFolder/Output").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathOutput->setText(path);
    ui->lineEditPathOutput->setCursorPosition( path.size() );
    m_strOutputFolder = path;
  }
  path = s.value("CurrentFolder/Calibration").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathCalibrationFile->setText(path);
    ui->lineEditPathCalibrationFile->setCursorPosition( path.size() );
    m_strCalibrationFile = path;
  }

  ui->widgetImageView->SetNumberOfExpectedPoints(4);

  m_wndProgress = new ProgressWindow(this);
  m_wndProgress->hide();

  qRegisterMetaType<QProcess::ExitStatus>("QProcess::ExitStatus");
  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()), Qt::QueuedConnection);
  connect(m_proc, SIGNAL(errorOccurred(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));

  connect(&m_fileWatcher, SIGNAL(directoryChanged(QString)), SLOT(OnFileChanged(QString)), Qt::QueuedConnection);
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("CurrentFolder/Input", m_strInputFolder);
  s.setValue("CurrentFolder/Output", m_strOutputFolder);
  s.setValue("CurrentFolder/FinalOutput", m_strFinalOutputFolder);
  s.setValue("CurrentFolder/Calibration", m_strCalibrationFile);
  s.setValue("MainWindow/Geometry", geometry());
  s.setValue("MainWindow/SplitterState",  ui->splitter->saveState());

  delete ui;
}

void MainWindow::resizeEvent(QResizeEvent* e)
{
  RepositionProgressWindow();
  QMainWindow::resizeEvent(e);
}

void MainWindow::moveEvent(QMoveEvent* e)
{
  RepositionProgressWindow();
  QMainWindow::moveEvent(e);
}


void MainWindow::RepositionProgressWindow()
{
  QRect rc = m_wndProgress->rect();
  rc.moveCenter(rect().center());
  rc.moveTopLeft(rc.topLeft()+this->pos());
  m_wndProgress->setGeometry(rc);
}

void MainWindow::OnButtonInputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->lineEditPathInput->text());
  if (!path.isEmpty())
  {
    ui->lineEditPathInput->setText(path);
    ui->lineEditPathInput->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonOutputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->lineEditPathOutput->text());
  if (!path.isEmpty())
  {
    ui->lineEditPathOutput->setText(path);
    ui->lineEditPathOutput->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonOutputFolderCC()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->labelPathCC->text());
  if (!path.isEmpty())
  {
    ui->labelPathCC->setText(path.trimmed());
    m_strFinalOutputFolder = path.trimmed();
    QFileInfoList list = QDir(m_strFinalOutputFolder).entryInfoList(QDir::Files, QDir::Name);
    bool bRun = true;
    if (list.size() == m_listInputFiles.size())
    {
      QMessageBox msgbox(this);
      msgbox.setIcon(QMessageBox::Question);
      QAbstractButton* yesBtn = msgbox.addButton("Overwrite", QMessageBox::YesRole);
      msgbox.addButton("Cancel", QMessageBox::NoRole);
      msgbox.setText("Existing files found in the folder. You can choose to overwrite them. Or cancel and select a different folder.");
      msgbox.setWindowTitle("Existing Files");
      msgbox.exec();
      if (msgbox.clickedButton() != yesBtn)
        bRun = false;
    }

    if (bRun)
    {
      QStringList cmd;
      cmd << m_strPythonCmd << m_strPyScriptMaskToCC
          << "--in_dir" << QString("\"%1\"").arg(m_strMaskFolder)
          << "--out_dir" << QString("\"%1\"").arg(m_sTempDir);
      m_proc->start(cmd.join(" "));
      m_proc->setProperty("task", "mask_to_cc");
      ui->pushButtonCreateMask->setEnabled(false);
      ShowProgressWindow(true);
    }
  }
}

void MainWindow::ShowProgressWindow(bool bShow)
{
  m_wndProgress->setVisible(bShow);
  if (bShow)
    RepositionProgressWindow();
}

void MainWindow::OnButtonCalibrationFile()
{
  QString path = QFileDialog::getOpenFileName(this, "Select Calibration File", ui->lineEditPathCalibrationFile->text());
  if (!path.isEmpty())
  {
    ui->lineEditPathCalibrationFile->setText(path);
    ui->lineEditPathCalibrationFile->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonLoadMask()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->lineEditPathOutput->text());
  if (!path.isEmpty())
  {
    ui->labelMaskFolder->setText(path);
    m_strMaskFolder = path;
    m_listMaskFiles = QDir(path).entryInfoList(QDir::Files, QDir::Name);
    m_listInputFiles = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
    if (!m_listInputFiles.isEmpty())
    {
      bool bRunNNUnet = true;
      if (m_listInputFiles.size() == m_listMaskFiles.size())
      {
        QMessageBox msgbox(this);
        msgbox.setIcon(QMessageBox::Question);
        QAbstractButton* yesBtn = msgbox.addButton("Re-generate and Overwrite", QMessageBox::YesRole);
        msgbox.addButton("Load Existing Masks", QMessageBox::NoRole);
        msgbox.setText("Existing masks found in the folder. You can choose to load them directly, or re-generate with nnunet. Re-generating will overwrite the existing files.");
        msgbox.setWindowTitle("Existing Masks");
        msgbox.exec();
        if (msgbox.clickedButton() != yesBtn)
          bRunNNUnet = false;
      }
      if (bRunNNUnet)
      {
        ClearFolder(path);
        QStringList cmd;
        QDir dir(m_sTempDir);
        QString temp_in = QString("tmp_in_%1").arg(QDateTime::currentMSecsSinceEpoch());
        QString temp_out = QString("tmp_out_%1").arg(QDateTime::currentMSecsSinceEpoch());
        dir.mkdir(temp_in);
        dir.mkdir(temp_out);
        QString temp_out_dir = m_sTempDir+"/"+temp_out;
        QStringList watch_list;
        m_mapInputImageSize.clear();
        for (int i = 0; i < m_listInputFiles.size(); i++)
        {
          QFile::copy(m_listInputFiles[i].absoluteFilePath(),
                      QFileInfo(QString("%1/%2").arg(m_sTempDir).arg(temp_in), m_listInputFiles[i].fileName()).absoluteFilePath());
          watch_list << QFileInfo(temp_out_dir, m_listInputFiles[i].completeBaseName()+".npz").absoluteFilePath();
          m_mapInputImageSize[watch_list.last()] = QImage(m_listInputFiles[i].absoluteFilePath()).size();
        }
        cmd << m_strPythonCmd << QString("\"%1/process_directory_no_upsample.py\"").arg(m_strNNUnetScriptFolder)
            << "--in_dir" << QString("\"%1\"").arg(m_sTempDir+"/"+temp_in)
            << "--out_dir" << QString("\"%1\"").arg(temp_out_dir)
            << "--model_path" << QString("\"%1\"").arg(m_strNNUnetModelFolder);
        m_proc->start(cmd.join(" "));
        m_proc->setProperty("task", "nnunet");
        m_proc->setProperty("output_folder", path);
        m_proc->setProperty("temp_output_folder", temp_out_dir);
        ui->pageSegEdit->setEnabled(false);
        m_elapsedTimer.start();
        m_fileWatcher.addPath(temp_out_dir);
        m_listQueuedFiles = watch_list;
        m_nIndex = 0;
        ShowProgressWindow(true);
      }
      else
      {
        LoadImage(0);
        ui->widgetSegCtrls->setEnabled(true);
        UpdateIndex();
      }
    }
  }
}

void MainWindow::OnButtonContinue()
{
  m_strInputFolder = ui->lineEditPathInput->text().trimmed();
  m_strOutputFolder = ui->lineEditPathOutput->text().trimmed();
  m_strCalibrationFile = ui->lineEditPathCalibrationFile->text().trimmed();

  if (m_strInputFolder.isEmpty() || !QDir(m_strInputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Input directory is not set or does not exist.");
    return;
  }
  if (m_strOutputFolder.isEmpty() || !QDir(m_strOutputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Output directory for corrected images is not set or does not exist.");
    return;
  }

  if (!m_strCalibrationFile.isEmpty() && !QFile(m_strOutputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Calibration file does not exist.");
    return;
  }
  m_listInputFiles = QDir(m_strInputFolder).entryInfoList(QDir::Files, QDir::Name);

  ui->stackedWidget->setCurrentWidget(ui->pageCorrection);
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_POINT);
  ui->widgetImageView->SetNumberOfExpectedPoints(4);
  ui->widgetImageView->SetMaskOpacity(0.7);

  m_bCalibratiedMode = !m_strCalibrationFile.isEmpty();
  ui->labelTitle->setText(m_bCalibratiedMode?"Calibrated Mode":"Retrospective Mode");

  if (!m_listInputFiles.isEmpty())
    LoadImage(m_nIndex);
  UpdateIndex();

  ui->pushButtonNext->setEnabled(false);

  QFileInfoList flist = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
  if (flist.size() == m_listInputFiles.size())
  {
    ui->pushButtonSegmentation->setEnabled(true);
    QMessageBox::information(this->ui->widgetControl, "Existing Corrected Images",
                             "Existing corrected images found in the folder. You can skip this ection by clicking 'Go to Segmentation' button directly.");
  }
}

void MainWindow::OnTogglePointMode(bool b)
{
  if (!b)
    return;

  bool b2Points = (sender() == ui->radioButton2Points);
  ui->labelRectWidth->setText(b2Points?"Ruler width (mm)":"Rectangle width (mm)");
  ui->labelRectHeight->setVisible(!b2Points);
  ui->lineEditRectHeight->setVisible(!b2Points);
  ui->widgetImageView->SetNumberOfExpectedPoints(b2Points?2:(sender() == ui->radioButton3Points?3:4));
}

void MainWindow::SetupScriptPath()
{  
  m_strNNUnetScriptFolder = QProcessEnvironment::systemEnvironment().value("NNUNET_SCRIPT_DIR");
  m_strNNUnetModelFolder = QProcessEnvironment::systemEnvironment().value("NNUNET_MODEL_DIR");
  if (m_strNNUnetModelFolder.isEmpty())
    qDebug() << "NNUNET_MODEL_DIR is not set!";
  if (m_strNNUnetScriptFolder.isEmpty())
    qDebug() << "NNUNET_SCRIPT_DIR is not set!";

  // copy resource files
  static QTemporaryDir dir;
  m_sTempDir = dir.path();
  // for retrospective correction
  {
    QFile file(QString(":/")+SCRIPT_RETROSPECTIVE);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QString str = in.readAll();
    str.replace("./resources/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
    str.replace("./resources/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
    file.close();
    file.flush();
    file.setFileName(QFileInfo(dir.path(),SCRIPT_RETROSPECTIVE).absoluteFilePath());
    if (file.open(QFile::ReadWrite))
    {
      QTextStream out(&file);
      out << str;
    }
    file.close();
    file.flush();
    QFile::copy(":/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
    QFile::copy(":/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
    m_strPyScriptRetrospective = file.fileName();

    if (!QFile::exists(m_strPyScriptRetrospective))
    {
      QMessageBox::critical(this, "Error", tr("Could not locate %1 script").arg(SCRIPT_RETROSPECTIVE));
      QTimer::singleShot(0, qApp, SLOT(quit()));
    }
  }
  // for calibration
  {
    QFile file(QString(":/")+SCRIPT_FIDUCIALS_CORRECTION);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QString str = in.readAll();
    str.replace("./resources/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
    str.replace("./resources/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
    file.close();
    file.flush();
    file.setFileName(QFileInfo(dir.path(),SCRIPT_FIDUCIALS_CORRECTION).absoluteFilePath());
    if (file.open(QFile::ReadWrite))
    {
      QTextStream out(&file);
      out << str;
    }
    file.close();
    file.flush();
    m_strPyScriptFiducialsCorrection = file.fileName();
    QFile::copy(":/registration.py", QFileInfo(dir.path(),"registration.py").absoluteFilePath());
    m_strPyScriptFiducialsDetection = QFileInfo(dir.path(),"func_fiducials_detection.py").absoluteFilePath();
    QFile::copy(":/func_fiducials_detection.py", m_strPyScriptFiducialsDetection);
    m_strPyScriptFiducialsCalibration = QFileInfo(dir.path(),"func_fiducials_calibration.py").absoluteFilePath();
    QFile::copy(":/func_fiducials_calibration.py", m_strPyScriptFiducialsCalibration);

    if (!QFile::exists(m_strPyScriptFiducialsCorrection))
    {
      QMessageBox::critical(this, "Error", tr("Could not locate %1 script").arg(SCRIPT_RETROSPECTIVE));
      QTimer::singleShot(0, qApp, SLOT(quit()));
    }
  }
  // for CC
  {
    m_listStockColors << QColor(255,100,100) << QColor(255,255,100) << QColor(100,255,100)
                      << QColor(110,245,255) << QColor(75,100,255) << QColor(255,128,0)
                      << QColor(100,150,170) << QColor(120,60,220);
    m_strPyScriptMaskToCC = QFileInfo(dir.path(),"func_mask_to_cc.py").absoluteFilePath();
    QFile::copy(":/func_mask_to_cc.py", m_strPyScriptMaskToCC);
  }
}

void MainWindow::UpdateIndex()
{
  QLabel* label = ui->labelIndex;
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    QFileInfoList flist = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
    ui->pushButtonPrevious->setEnabled(m_nIndex > 0);
    ui->pushButtonNext->setEnabled(m_nIndex < flist.size());
    ui->pushButtonSegmentation->setEnabled(flist.size() == m_listInputFiles.size());
  }
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    label = ui->labelIndexSeg;
    m_listMaskFiles = QDir(m_strMaskFolder).entryInfoList(QDir::Files, QDir::Name);
    ui->pushButtonPreviousSeg->setEnabled(m_nIndex > 0);
    ui->pushButtonNextSeg->setEnabled(m_nIndex < m_listMaskFiles.size()-1);
    ui->pushButtonCC->setEnabled(m_listInputFiles.size() == m_listMaskFiles.size());
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    label = ui->labelIndexCC;
    ui->pushButtonPreviousCC->setEnabled(m_nIndex > 0);
    ui->pushButtonNextCC->setEnabled(m_nIndex < m_listRegionData.size());
  }
  label->setText(tr("%1 / %2").arg(m_nIndex+1).arg(m_listInputFiles.size()));
}

void MainWindow::LoadImage(int n)
{
  QList<QPoint> pts;
  QList<RECT_REGION> rects;
  QString fn = m_listInputFiles[n].absoluteFilePath(), mask_fn;
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    if (m_listPointData.size() > n)
      pts = m_listPointData[n];
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    m_listMaskFiles = QDir(m_strMaskFolder).entryInfoList(QDir::Files, QDir::Name);
    mask_fn = m_listMaskFiles[n].absoluteFilePath();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    mask_fn = m_listMaskFiles[n].absoluteFilePath();
    if (m_listRegionData.size() > n)
      rects = m_listRegionData[n];
  }
  ui->widgetImageView->LoadImage(fn, mask_fn, pts, rects);
  if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    if (pts.isEmpty())
    {
      QStringList cmd;
      cmd << m_strPythonCmd << m_strPyScriptFiducialsDetection
          << "--in_image" << QString("\"%1\"").arg(fn)
          << "--calibration_file" << QString("\"%1\"").arg(m_strCalibrationFile)
          << "--out_file" << m_strPyScriptFiducialsDetection + ".txt";
      m_proc->start(cmd.join(" "));
      m_proc->setProperty("task", "fiducials_detection");
      m_proc->setProperty("output_file", cmd.last());
      ui->pageCorrection->setEnabled(false);
    }
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    QString fn_only = m_listMaskFiles[n].fileName();
    fn_only.replace(QString(".")+QFileInfo(mask_fn).suffix(), ".npz");
    if (!m_maskProcessor.Load(QFileInfo(m_sTempDir, fn_only).absoluteFilePath()))
    {
      QMessageBox::warning(this, "Error", "Failed to load processed npy files");
      return;
    }

    if (!rects.isEmpty())
    {
      QList<QPoint> pts;
      for (int i = 0; i < rects.size(); i++)
        pts << rects[i].first << rects[i].second;
      m_maskProcessor.LoadSelections(pts);
      ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
    }
  }
}

QList<QPoint> MainWindow::GetCalibrationPointsList(const QVariantMap& info)
{
  QList<QPoint> pts;
  QPoint pt;
  for (int i = 0; i < 4; i++)
  {
    QString id_str = QString("corner%1").arg(i);
    pt.setX(info[id_str].toList().first().toDouble());
    pt.setY(info[id_str].toList().last().toDouble());
    pts << pt;
  }
  return pts;
}

void MainWindow::OnButtonNext()
{
  if (m_nIndex < m_listInputFiles.size()-1)
  {
    if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
    {
      ui->widgetImageView->SaveMaskIfEdited();
    }
    LoadImage(++m_nIndex);
    UpdateIndex();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    ui->pushButtonAllDone->setEnabled(true);
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCorrection)
  {
    ui->pushButtonSegmentation->setEnabled(true);
  }
}

void MainWindow::OnButtonPrevious()
{
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    ui->widgetImageView->SaveMaskIfEdited();
  }
  LoadImage(--m_nIndex);
  UpdateIndex();
}

void MainWindow::OnButtonProcess()
{
  bool bOK;
  double dWidth = ui->lineEditRectWidth->text().trimmed().toDouble(&bOK);
  double dHeight = 1;
  int nExpectedPoints = 2;
  if (ui->radioButton4Points->isChecked() || ui->radioButton3Points->isChecked())
  {
    nExpectedPoints = (ui->radioButton4Points->isChecked()?4:3);
    dHeight = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
      dHeight = -1;
  }
  if (false)
  {
    double w = 1, h = 1;
    bool bOK;
    w = ui->lineEditRectWidth->text().trimmed().toDouble(&bOK);
    if (bOK)
      h = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
    {
      QMessageBox::warning(this, "Error", "Please enter valid width and height");
      return;
    }
    QStringList cmd;
    QStringList list;
    QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
    foreach (QPoint pt, pts)
      list << QString::number(pt.x()) << QString::number(pt.y());
    QString fn = QFileInfo(ui->widgetImageView->GetFilename()).fileName();
    QString out_fn = QFileInfo(m_sTempDir, fn.left(fn.lastIndexOf('.'))+"_corrected.npz").absoluteFilePath();
    cmd << m_strPythonCmd << m_strPyScriptFiducialsCalibration
        << "--in_img" << QString("\"%1\"").arg(ui->widgetImageView->GetFilename())
        << "--points" << list.join(" ")
        << "--width" << QString::number(w)
        << "--height" << QString::number(h)
        << "--out_file" << QString("\"%1\"").arg(out_fn);
    m_proc->start(cmd.join(" "));
    m_proc->setProperty("task", "fiducials_calibration");
    m_proc->setProperty("output_file", out_fn);
  }
  else if (ui->widgetImageView->GetEditedPoints().size() < nExpectedPoints)
    QMessageBox::warning(this, "", "Please click on the image to add another point");
  else if (!bOK || dWidth <= 0 || dHeight <= 0)
    QMessageBox::warning(this, "", "Please enter a valid value for rectangle length");
  else
  {
    // execute script
    QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
    if (m_nIndex < m_listPointData.size())
      m_listPointData[m_nIndex] = pts;
    else
      m_listPointData << pts;

    QStringList strList;
    foreach (QPoint pt, pts)
      strList << QString::number(pt.x()) << QString::number(pt.y());
    QStringList cmd;
    cmd << m_strPythonCmd << m_strPyScriptRetrospective
        << "--in_img" << QString("\"%1\"").arg(ui->widgetImageView->GetFilename())
        << "--points" << strList.join(" ")
        << "--width" << QString::number(dWidth)
        << "--height" << QString::number(dHeight)
        << "--out_dir" << QString("\"%1\"").arg(m_strOutputFolder);
    m_proc->start(cmd.join(" "));
    m_proc->setProperty("task", "retrospective_correction");
  }
}

void MainWindow::OnProcessError(QProcess::ProcessError er)
{
  QString str = "Unknown error occurred during process.";
  switch (er)
  {
  case QProcess::FailedToStart:
    str = "Failed to start python script";
    break;
  case QProcess::Crashed:
    str = "Process crashed";
    break;
  case QProcess::Timedout:
    str = "Process timed out";
    break;
  case QProcess::ReadError:
    str = "Failed to read";
    break;
  case QProcess::WriteError:
    str = "Failed to write";
    break;
  default:
    break;
  }

  ShowProgressWindow(false);
  qDebug() << str;
}

void MainWindow::OnProcessStarted()
{
  ui->pushButtonProcess->setEnabled(false);
  ui->widgetImageView->ShowMessage("Processing...");
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->HideMessage();
  if (m_proc->exitStatus() == QProcess::NormalExit)
    ui->pushButtonNext->setEnabled(true);
  ui->pushButtonProcess->setEnabled(true);
  if (m_nIndex >= m_listInputFiles.size()-1)
    ui->pushButtonSegmentation->setEnabled(true);

  ShowProgressWindow(false);

  QString task;
  if (sender())
    task = sender()->property("task").toString();
  if (task == "fiducials_correction")
  {
    QStringList cmd;
    cmd << m_strPythonCmd << m_strPyScriptFiducialsDetection
        << "--in_image" << ui->widgetImageView->GetFilename()
        << "--calibration_file" << QString("\"%1\"").arg(ui->lineEditPathCalibrationFile->text().trimmed())
        << "--out_file" << m_strPyScriptFiducialsDetection + ".txt";
    m_proc->start(cmd.join(" "));
    m_proc->setProperty("task", "fiducials_detection");
    m_proc->setProperty("output_file", cmd.last());
    ui->widgetImageView->SetEditMode(WidgetImageView::EM_CALIBRATION);
  }
  else if (task == "fiducials_detection")
  {
    QFile file(sender()->property("output_file").toString());
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QString str = in.readAll();
    QStringList list = str.split("\n");
    if (list.size() >= 6)
    {
      QVariantMap info;
      info["width"] = list[0].split(":").last().trimmed().toDouble();
      info["height"] = list[1].split(":").last().trimmed().toDouble();
      info["corner0"] = list[2].split(":").last().trimmed().split(",");
      info["corner1"] = list[3].split(":").last().trimmed().split(",");
      info["corner2"] = list[4].split(":").last().trimmed().split(",");
      info["corner3"] = list[5].split(":").last().trimmed().split(",");

      ui->lineEditRectWidth->setText(info["width"].toString());
      ui->lineEditRectHeight->setText(info["height"].toString());
      QList<QPoint> pts = GetCalibrationPointsList(info);
      ui->widgetImageView->SetEditedPoints(pts);
    }
    ui->pageCorrection->setEnabled(true);
    UpdateIndex();
  }
  else if (task == "fiducials_calibration")
  {
    QString calibration_file = sender()->property("output_file").toString();
    //    qDebug() << calibration_file;
    QDir dir(m_sTempDir);
    QString fn = ui->widgetImageView->GetFilename();
    QString sub_dir = QFileInfo(fn).completeBaseName();
    dir.mkdir(sub_dir);
    dir.cd(sub_dir);
    QFile::copy(fn, QFileInfo(dir.path(), QFileInfo(fn).fileName()).absoluteFilePath());
    QStringList cmd;
    cmd << m_strPythonCmd << m_strPyScriptFiducialsCorrection
        << "--in_dir" << QString("\"%1\"").arg(dir.path())
        << "--calibration_file" << QString("\"%1\"").arg(calibration_file)
        << "--out_dir" << QString("\"%1\"").arg(m_strOutputFolder);
    m_proc->start(cmd.join(" "));
    m_proc->setProperty("task", "fiducials_correction_single");
  }
  else if (task == "mask_to_cc")
  {
    if (!m_listInputFiles.isEmpty() && !m_listMaskFiles.isEmpty())
      LoadImage(0);

    ui->widgetImageView->HideMessage();
    ui->pushButtonCreateMask->setEnabled(true);
    UpdateIndex();
  }
  else if (task == "nnunet")
  {
    qDebug() << "nnUNet elapsed time in secs: " << m_elapsedTimer.elapsed()/1000;
    ui->pageSegEdit->setEnabled(true);
    UpdateIndex();
    OnFileChanged(m_proc->property("temp_output_folder").toString());
  }
}

void MainWindow::OnProcessOutputMessage()
{
  qDebug() << m_proc->readAllStandardOutput();
}

void MainWindow::OnProcessErrorMessage()
{
  QString str = m_proc->readAllStandardError().trimmed();
  if (!str.isEmpty())
  {
    qDebug() << str;
  }
}

void MainWindow::OnButtonClear()
{
  ui->widgetImageView->ClearEdits();
  m_maskProcessor.ClearBuffer();
  UpdateIndex();
}

void MainWindow::OnButtonProceedToSeg()
{
  ui->stackedWidget->setCurrentWidget(ui->pageSegEdit);
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_EDIT_MASK);
  m_nIndex = 0;
}

void MainWindow::OnButtonProceedToCC()
{
  ui->widgetImageView->SaveMaskIfEdited();

  ui->stackedWidget->setCurrentWidget(ui->pageCC);
  m_nIndex = 0;
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_REGION);
  connect(ui->widgetImageView, SIGNAL(LastRegionEdited(int)),
          SLOT(OnLastRegionEdited(int)), Qt::QueuedConnection);

  m_listInputFiles = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
}

void MainWindow::OnSliderSegOpacity(int n)
{
  ui->widgetImageView->SetMaskOpacity(n/100.0);
}

void MainWindow::OnLastRegionEdited(int n)
{
  RECT_REGION rc = ui->widgetImageView->GetEditedRegions().last();
  if (m_maskProcessor.ProcessSelection(rc.first, rc.second, n+1))
  {
    ui->widgetImageView->SetOverlay(m_maskProcessor.GetMaskImage(m_listStockColors));
  }
  else
  {
    QMessageBox::warning(this, "Error", "Did not find any new slice");
  }
}

void MainWindow::OnButtonCreateMask()
{
  QList<RECT_REGION> list = ui->widgetImageView->GetEditedRegions();
  if (m_nIndex < m_listRegionData.size())
    m_listRegionData[m_nIndex] = list;
  else
    m_listRegionData << list;

  QString fn = QFileInfo(ui->widgetImageView->GetFilename()).fileName();
  fn.replace(QString(".")+QFileInfo(fn).suffix(), "_mask.npy", Qt::CaseInsensitive);
  m_maskProcessor.SaveToNpy(QFileInfo(m_strFinalOutputFolder, fn).absoluteFilePath());
  ui->pushButtonNextCC->setEnabled(true);
  ui->widgetImageView->ShowMessage("Mask created and saved");
}

void MainWindow::OnToggleMask()
{
  if (ui->widgetImageView->GetMaskOpacity() == 0)
    ui->widgetImageView->SetMaskOpacity(ui->horizontalSliderSegOpacity->value()/100.0);
  else
    ui->widgetImageView->SetMaskOpacity(0);
}

void MainWindow::OnFileChanged(const QString& path)
{
  static QVariantMap last_size_info = QVariantMap();
  QFileInfoList flist = QDir(path).entryInfoList(QStringList("*.npz"), QDir::Files, QDir::Name);
  while (!flist.isEmpty())
  {
    if (!flist.isEmpty() && !m_listQueuedFiles.isEmpty())
    {
      QString fn = flist[0].absoluteFilePath();
      if (fn == m_listQueuedFiles.first())
      {
        if (last_size_info[fn].toInt() == flist[0].size() && flist[0].size() > 0)
        {
          QImage image = NpyToImage(fn, m_mapInputImageSize[fn].toSize());
          QString fn_png = QFileInfo(m_strMaskFolder, QFileInfo(fn).completeBaseName()+".png").absoluteFilePath();
          image.save(fn_png);
          m_listQueuedFiles.removeFirst();
          QFile::remove(fn);
          ui->pageSegEdit->setEnabled(true);
          ui->widgetSegCtrls->setEnabled(true);
          if (ui->widgetImageView->GetMaskFilename().isEmpty())
          {
            ShowProgressWindow(false);
            LoadImage(0);
            QMessageBox::information(this, "Edit Mask", "You can start editing the current mask while waiting for the rest to be generated.");
          }
          UpdateIndex();
          flist = QDir(path).entryInfoList(QStringList("*.npz"), QDir::Files, QDir::Name);
        }
        else
        {
          last_size_info[fn] = flist[0].size();
          break;
        }
      }
      else
        break;
    }
    else
      break;
  }
}

QImage MainWindow::NpyToImage(const QString& npy_in, const QSize& sz_in)
{
  QImage image;
  if (!QFile::exists(npy_in))
  {
    qDebug() << "File does not exist:" << npy_in;
    return image;
  }
  cnpy::NpyArray ar = cnpy::npz_load(qPrintable(npy_in), "probabilities");
  if (ar.shape.size() != 4)
  {
    qDebug() << "Could not load numpy file " << npy_in;
    return image;
  }

  float* ptr = ar.data<float>();
  int w = ar.shape[3];
  int h = ar.shape[2];
  image = QImage(QSize(w, h), QImage::Format_ARGB32);
  image.fill(QColor(0,0,0,255));

  for (int y = 0; y < h; y++)
  {
    QRgb* p = (QRgb*)image.scanLine(y);
    for (int x = 0; x < w; x++)
    {
      if (ptr[w*h+y*w+x] >= 0.5)
      {
        p[x] = qRgba(255,255,255,255);
      }
    }
  }

  QSize sz = sz_in;
  if (!sz.isValid())
    sz = QImage(m_listInputFiles.first().absoluteFilePath()).size();
  image = image.scaled(sz, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
  w = image.width();
  h = image.height();
  for (int y = 0; y < h; y++)
  {
    QRgb* p = (QRgb*)image.scanLine(y);
    for (int x = 0; x < w; x++)
    {
      if (qRed(p[x]) >= 128)
        p[x] = qRgba(255,255,255,255);
      else
        p[x] = qRgba(0,0,0,255);
    }
  }

  return image;
}

void MainWindow::ClearFolder(const QString& path)
{
  QFileInfoList flist = QDir(path).entryInfoList(QDir::Files, QDir::Name);
  foreach (QFileInfo fi, flist)
    QFile::remove(fi.absoluteFilePath());
}
