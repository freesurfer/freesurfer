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

#define SCRIPT_RETROSPECTIVE "func_retrospective_correction.py"
#define SCRIPT_FIDUCIALS_CORRECTION "func_fiducials_correction.py"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);
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
  path = s.value("CurrentFolder/FinalOutput").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathFinalOutput->setText(path);
    ui->lineEditPathFinalOutput->setCursorPosition( path.size() );
    m_strFinalOutputFolder = path;
  }
  path = s.value("CurrentFolder/Calibration").toString();
  if (!path.isEmpty())
  {
    ui->lineEditPathCalibrationFile->setText(path);
    ui->lineEditPathCalibrationFile->setCursorPosition( path.size() );
    m_strCalibrationFile = path;
  }

  ui->widgetImageView->SetNumberOfExpectedPoints(4);

  qRegisterMetaType<QProcess::ExitStatus>("QProcess::ExitStatus");
  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()), Qt::QueuedConnection);
  connect(m_proc, SIGNAL(errorOccurred(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));
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

void MainWindow::OnButtonFinalOutputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->lineEditPathFinalOutput->text());
  if (!path.isEmpty())
  {
    ui->lineEditPathFinalOutput->setText(path);
    ui->lineEditPathFinalOutput->setCursorPosition(path.size());
  }
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
  QString path = QFileDialog::getExistingDirectory(this, "Select Folder", ui->lineEditPathFinalOutput->text());
  if (!path.isEmpty())
  {
    ui->labelMaskFolder->setText(path);
    m_strMaskFolder = path;
    m_listMaskFiles = QDir(path).entryInfoList(QDir::Files, QDir::Name);
    m_listInputFiles = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
    if (!m_listInputFiles.isEmpty() && !m_listMaskFiles.isEmpty())
    {
      if (m_listInputFiles.size() != m_listMaskFiles.size())
      {
        qDebug() << "File counts do not match";
      }
      LoadImage(0);
      ui->widgetSegCtrls->setEnabled(true);
      UpdateIndex();
    }
  }
}

void MainWindow::OnButtonContinue()
{
  m_strInputFolder = ui->lineEditPathInput->text().trimmed();
  m_strOutputFolder = ui->lineEditPathOutput->text().trimmed();
  m_strFinalOutputFolder = ui->lineEditPathFinalOutput->text().trimmed();
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
  if (m_strFinalOutputFolder.isEmpty() || !QDir(m_strFinalOutputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Output directory for connected components is not set or does not exist.");
    return;
  }
  if (!m_strCalibrationFile.isEmpty() && !QFile(m_strOutputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Calibration file does not exist.");
    return;
  }
  m_listInputFiles = QDir(m_strInputFolder).entryInfoList(QDir::Files, QDir::Name);

  ui->stackedWidget->setCurrentWidget(ui->pageCorrection);

  m_bCalibratiedMode = !m_strCalibrationFile.isEmpty();
  ui->labelTitle->setText(m_bCalibratiedMode?"Calibrated Mode":"Retrospective Mode");
  ui->widgetPointMode->setVisible(!m_bCalibratiedMode);
  if (m_bCalibratiedMode)
  {
    QStringList cmd;
    cmd << m_strPythonCmd << m_strPyScriptFiducialsCorrection
        << "--in_dir" << QString("\"%1\"").arg(m_strInputFolder)
        << "--calibration_file" << QString("\"%1\"").arg(ui->lineEditPathCalibrationFile->text().trimmed())
        << "--out_dir" << QString("\"%1\"").arg(m_strOutputFolder);
    m_proc->start(cmd.join(" "));
    m_proc->setProperty("task", "fiducials_correction");
    ui->pageCorrection->setEnabled(false);
  }

  UpdateIndex();

  if (!m_listInputFiles.isEmpty())
    LoadImage(m_nIndex);
}

void MainWindow::OnTogglePointMode(bool b)
{
  if (!b)
    return;

  bool b4Points = (sender() == ui->radioButton4Points);
  ui->labelRectWidth->setText(b4Points?"Rectangle width (mm)":"Ruler width (mm)");
  ui->labelRectHeight->setVisible(b4Points);
  ui->lineEditRectHeight->setVisible(b4Points);
  ui->widgetImageView->SetNumberOfExpectedPoints(b4Points?4:2);
}

void MainWindow::SetupScriptPath()
{
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
    ui->pushButtonPrevious->setEnabled(m_nIndex > 0);
    if (!m_bCalibratiedMode)
      ui->pushButtonNext->setEnabled(m_nIndex < m_listPointData.size());
  }
  if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    label = ui->labelIndexSeg;
    ui->pushButtonPreviousSeg->setEnabled(m_nIndex > 0);
    ui->pushButtonNextSeg->setEnabled(m_nIndex < m_listInputFiles.size()-1);
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    label = ui->labelIndexCC;
    ui->pushButtonPreviousSeg->setEnabled(m_nIndex > 0);
    ui->pushButtonNext->setEnabled(m_nIndex < m_listRegionData.size());
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
    if (!m_bCalibratiedMode)
    {
      if (m_listPointData.size() > n)
        pts = m_listPointData[n];
    }
    else
    {
      QVariantMap info;
      if (m_mapCalibrationInfo.contains(fn))
        info = m_mapCalibrationInfo[fn].toMap();
      else if (m_mapCalibrationInfo.contains("general"))
        info = m_mapCalibrationInfo["general"].toMap();
      if (!info.isEmpty())
      {
        pts = GetCalibrationPointsList(info);
      }
    }
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageSegEdit)
  {
    mask_fn = m_listMaskFiles[n].absoluteFilePath();
  }
  else if (ui->stackedWidget->currentWidget() == ui->pageCC)
  {
    mask_fn = m_listMaskFiles[n].absoluteFilePath();
    if (m_listRegionData.size() > n)
      rects = m_listRegionData[n];
  }
  ui->widgetImageView->LoadImage(fn, mask_fn, pts, rects);
  if (ui->stackedWidget->currentWidget() == ui->pageCC)
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
      ui->pushButtonNext->setEnabled(true);
    }
  }
}

QList<QPoint> MainWindow::GetCalibrationPointsList(const QVariantMap& info)
{
  QList<QPoint> pts;
  QPoint pt;
  for (int i = 0; i < 4; i++)
  {
    QString id_str = QString("corner%1").arg(i), id_str_2 = QString("end_point%1").arg(i);
    pt.setX(info[id_str].toList().first().toDouble());
    pt.setY(info[id_str].toList().last().toDouble());
    pts << pt;
    if (info.contains(id_str_2))
    {
      pt.setX(info[id_str_2].toList().first().toDouble());
      pt.setY(info[id_str_2].toList().last().toDouble());
    }
    else
    {
      pt.setX(pt.x()+2);
    }
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
  if (ui->radioButton4Points->isChecked() || m_bCalibratiedMode)
  {
    nExpectedPoints = 4;
    dHeight = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
      dHeight = -1;
  }
  if (m_bCalibratiedMode)
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
      m_mapCalibrationInfo["general"] = info;
      QString fn = ui->widgetImageView->GetFilename();
      if (!m_mapCalibrationInfo.contains(fn))
      {
        QList<QPoint> pts = GetCalibrationPointsList(info);
        ui->widgetImageView->SetEditedPoints(pts);
      }
    }
    ui->pageCorrection->setEnabled(true);
  }
  else if (task == "fiducials_calibration")
  {
    QString calibration_file = sender()->property("output_file").toString();
    qDebug() << calibration_file;
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
  ui->pushButtonNext->setEnabled(false);
  ui->pushButtonNextCC->setEnabled(false);
  ui->pushButtonAllDone->setEnabled(false);
  m_maskProcessor.ClearBuffer();
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
  ui->widgetImageView->SetEditMode(WidgetImageView::EM_REGION);
  connect(ui->widgetImageView, SIGNAL(LastRegionEdited(int)),
                                      SLOT(OnLastRegionEdited(int)), Qt::QueuedConnection);

  m_listInputFiles = QDir(m_strOutputFolder).entryInfoList(QDir::Files, QDir::Name);
  QStringList cmd;
  cmd << m_strPythonCmd << m_strPyScriptMaskToCC
      << "--in_dir" << QString("\"%1\"").arg(m_strMaskFolder)
      << "--out_dir" << QString("\"%1\"").arg(m_sTempDir);
  m_proc->start(cmd.join(" "));
  m_proc->setProperty("task", "mask_to_cc");
  ui->pushButtonCreateMask->setEnabled(false);
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
}

