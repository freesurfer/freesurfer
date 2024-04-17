#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "CommonDef.h"
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <QScreen>
#include <QTemporaryDir>
#include <QTimer>

#define SCRIPT_RETROSPECTIVE "func_retrospective_correction.py"


MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);
  ui->stackedWidget->setCurrentWidget(ui->page1);

  m_strPythonCmd = PY_COMMAND;

  SetupScriptPath();

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QGuiApplication::primaryScreen()->geometry().contains(rc))
    setGeometry(rc);

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

  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()));
  connect(m_proc, SIGNAL(errorOccurred(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("CurrentFolder/Input", m_strInputFolder);
  s.setValue("CurrentFolder/Output", m_strOutputFolder);
  s.setValue("CurrentFolder/Calibration", m_strCalibrationFile);

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

void MainWindow::OnButtonCalibrationFile()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Calibration File", ui->lineEditPathCalibrationFile->text());
  if (!path.isEmpty())
  {
    ui->lineEditPathCalibrationFile->setText(path);
    ui->lineEditPathCalibrationFile->setCursorPosition(path.size());
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
    QMessageBox::warning(this, "Error", "Output directory is not set or does not exist.");
    return;
  }
  if (!m_strCalibrationFile.isEmpty() && !QFile(m_strOutputFolder).exists())
  {
    QMessageBox::warning(this, "Error", "Calibration file does not exist.");
    return;
  }
  m_listInputFiles = QDir(m_strInputFolder).entryInfoList(QDir::Files, QDir::Name);

  ui->stackedWidget->setCurrentWidget(ui->page2);

  bool bRetroMode = m_strCalibrationFile.isEmpty();
  ui->labelTitle->setText(bRetroMode?"Retrospective Mode":"Calibrated Mode");
  ui->widgetPointMode->setVisible(bRetroMode);

  UpdateIndex();

  if (!m_listInputFiles.isEmpty())
    LoadImage(m_nIndex);
}

void MainWindow::OnTogglePointMode(bool b)
{
  bool b4Points = (b && sender() == ui->radioButton4Points);
  ui->labelRectWidth->setText(b4Points?"Rectangle width (mm)":"Ruler width (mm)");
  ui->labelRectHeight->setVisible(b4Points);
  ui->lineEditRectHeight->setVisible(b4Points);
}

void MainWindow::SetupScriptPath()
{
  // copy resource files
  static QTemporaryDir dir;
  QFile file(QString(":/")+SCRIPT_RETROSPECTIVE);
  file.open(QFile::ReadOnly | QFile::Text);
  QTextStream in(&file);
  QString str = in.readAll();
  str.replace("./resources/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
  str.replace("./resources/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
  QFile file2(QFileInfo(dir.path(),SCRIPT_RETROSPECTIVE).absoluteFilePath());
  if (file2.open(QFile::ReadWrite))
  {
    QTextStream out(&file2);
    out << str;
  }
  file2.close();
  file2.flush();
  QFile::copy(":/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
  QFile::copy(":/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
  m_strPyScriptRetrospective = QFileInfo(dir.path(),SCRIPT_RETROSPECTIVE).absoluteFilePath();

  if (!QFile::exists(m_strPyScriptRetrospective))
  {
    QMessageBox::critical(this, "Error", tr("Could not locate %1 script").arg(SCRIPT_RETROSPECTIVE));
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }
}

void MainWindow::UpdateIndex()
{
  ui->labelIndex->setText(tr("%1 / %2").arg(m_nIndex+1).arg(m_listInputFiles.size()));
  if (ui->stackedWidget->currentWidget() == ui->page2)
  {
    ui->pushButtonPrevious->setEnabled(m_nIndex > 0);
    ui->pushButtonNext->setEnabled(m_nIndex < m_listData.size());
  }
}

void MainWindow::LoadImage(int n)
{
  QList<QPoint> pts;
  if (m_listData.size() > n)
    pts = m_listData[n];
  ui->widgetImageView->LoadImage(m_listInputFiles[n].absoluteFilePath(), "", pts);
}

void MainWindow::OnButtonNext()
{
  if (m_nIndex < m_listInputFiles.size()-1)
  {
    LoadImage(++m_nIndex);
    UpdateIndex();
  }
}

void MainWindow::OnButtonPrevious()
{
  LoadImage(--m_nIndex);
  UpdateIndex();
}

void MainWindow::OnButtonProcess()
{
  bool bRetro = m_strCalibrationFile.isEmpty();
  bool bOK;
  double dWidth = ui->lineEditRectWidth->text().trimmed().toDouble(&bOK);
  double dHeight = 1;
  int nExpectedPoints = 2;
  if (ui->radioButton4Points->isChecked())
  {
    nExpectedPoints = 4;
    dHeight = ui->lineEditRectHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
      dHeight = -1;
  }
  if (ui->widgetImageView->GetEditedPoints().size() < nExpectedPoints)
    QMessageBox::warning(this, "", "Please click on the image to add another point");
  else if (!bOK || dWidth <= 0 || dHeight <= 0)
    QMessageBox::warning(this, "", "Please enter a valid value for rectangle length");
  else
  {
    // execute script
    QList<QPoint> pts = ui->widgetImageView->GetEditedPoints();
    if (m_nIndex < m_listData.size())
      m_listData[m_nIndex] = pts;
    else
      m_listData << pts;

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

}
