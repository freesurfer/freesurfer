#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "CommonDef.h"
#include <QTemporaryDir>
#include <QSettings>
#include <QFileInfo>
#include <QMessageBox>
#include <QDesktopWidget>
#include <QTimer>
#include <QFileDialog>
#include <QDebug>

#define WND_TITLE "Correction With Fiducials"
#define SCRIPT_FILENAME "func_fiducials_correction.py"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  m_strPythonCmd = PY_COMMAND;

  // setup script
  SetupScriptPath();

  setWindowTitle(WND_TITLE);

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QApplication::desktop()->screenGeometry(this).contains(rc))
    setGeometry(rc);

  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(started()), SLOT(UpdateWidgets()), Qt::QueuedConnection);
  connect(m_proc, SIGNAL(finished(int,QProcess::ExitStatus)), SLOT(UpdateWidgets()));
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(finished(int,QProcess::ExitStatus)), SLOT(OnProcessFinished()));

  UpdateWidgets();
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("CurrentFolder/Input", ui->lineEditInputFolder->text().trimmed());
  s.setValue("CurrentFolder/Output", ui->lineEditOutputFolder->text().trimmed());
  s.setValue("CurrentFolder/CalibrationFile", ui->lineEditCalibrationFile->text().trimmed());

  s.setValue("MainWindow/Geometry", geometry());

  if (m_proc && m_proc->state() == QProcess::Running)
  {
    m_proc->waitForFinished();
  }

  delete ui;
}

void MainWindow::showEvent(QShowEvent* e)
{
  QSettings s;
  QString path = s.value("CurrentFolder/Input").toString();
  if (!path.isEmpty())
  {
    ui->lineEditInputFolder->setText(path);
    ui->lineEditInputFolder->setCursorPosition( path.size() );
  }
  path = s.value("CurrentFolder/Output").toString();
  if (!path.isEmpty())
  {
    ui->lineEditOutputFolder->setText(path);
    ui->lineEditOutputFolder->setCursorPosition( path.size() );
  }
  path = s.value("CurrentFolder/CalibrationFile").toString();
  if (!path.isEmpty())
  {
    ui->lineEditCalibrationFile->setText(path);
    ui->lineEditCalibrationFile->setCursorPosition( path.size() );
  }
}

void MainWindow::SetupScriptPath()
{
  // copy resource files
  static QTemporaryDir dir;
  QFile file(QString(":/")+SCRIPT_FILENAME);
  file.open(QFile::ReadOnly | QFile::Text);
  QTextStream in(&file);
  QString str = in.readAll();
  str.replace("./resources/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
  str.replace("./resources/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
  QFile file2(QFileInfo(dir.path(),SCRIPT_FILENAME).absoluteFilePath());
  if (file2.open(QFile::ReadWrite))
  {
    QTextStream out(&file2);
    out << str;
  }
  file2.close();
  file2.flush();
  QFile::copy(":/horizontal.png", QFileInfo(dir.path(),"horizontal.png").absoluteFilePath());
  QFile::copy(":/vertical.png", QFileInfo(dir.path(),"vertical.png").absoluteFilePath());
  m_strPyScriptPath = QFileInfo(dir.path(),SCRIPT_FILENAME).absoluteFilePath();
  QFile::copy(":/registration.py", QFileInfo(dir.path(),"registration.py").absoluteFilePath());

  if (!QFile::exists(m_strPyScriptPath))
  {
    QMessageBox::critical(this, "Error", tr("Could not locate %1 script").arg(SCRIPT_FILENAME));
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }
}


void MainWindow::UpdateWidgets()
{
  if (ui->lineEditCalibrationFile->text().trimmed().isEmpty() ||
      ui->lineEditInputFolder->text().trimmed().isEmpty() ||
      ui->lineEditOutputFolder->text().trimmed().isEmpty() ||
      m_proc->state() == QProcess::Running)
    ui->pushButtonGo->setEnabled(false);
  else
    ui->pushButtonGo->setEnabled(true);

  ui->labelMessage->setText(m_proc->state() == QProcess::Running?"Processing...":"");
}

void MainWindow::OnButtonInputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Input Folder", ui->lineEditInputFolder->text().trimmed());
  if (!path.isEmpty())
  {
    ui->lineEditInputFolder->setText(path);
    ui->lineEditInputFolder->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonOutputFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, "Select Output Folder", ui->lineEditOutputFolder->text().trimmed());
  if (!path.isEmpty())
  {
    ui->lineEditOutputFolder->setText(path);
    ui->lineEditOutputFolder->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonCalibrationFile()
{
  QString path = QFileDialog::getOpenFileName(this, "Select Calibration File", ui->lineEditCalibrationFile->text().trimmed());
  if (!path.isEmpty())
  {
    ui->lineEditCalibrationFile->setText(path);
    ui->lineEditCalibrationFile->setCursorPosition(path.size());
  }
}

void MainWindow::OnButtonProcess()
{
  QStringList cmd;
  cmd << m_strPythonCmd << m_strPyScriptPath
      << "--in_dir" << ui->lineEditInputFolder->text().trimmed()
      << "--calibration_file" << ui->lineEditCalibrationFile->text().trimmed()
      << "--out_dir" << ui->lineEditOutputFolder->text().trimmed();
//  qDebug() << cmd.join(" ");
  m_proc->start(cmd.join(" "));
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

void MainWindow::OnProcessFinished()
{
  QMessageBox dlg(QMessageBox::Information, "Correction", "Correction finished. Click Exit to quit.", QMessageBox::Ok, this);
  dlg.setButtonText(QMessageBox::Ok, "Exit");
  dlg.exec();
  QTimer::singleShot(0, qApp, SLOT(quit()));
}
