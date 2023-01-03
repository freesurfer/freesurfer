#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "DialogWelcome.h"
#include "DialogSelectFolder.h"
#include <QDebug>
#include <QTimer>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QSettings>
#include <QDesktopWidget>
#include "CommonDef.h"
#include <QTemporaryDir>

#define WND_TITLE "Retrospective Correction"
#define SCRIPT_FILENAME "func_retrospective_correction.py"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);

  m_strPythonCmd = PY_COMMAND;

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QApplication::desktop()->screenGeometry(this).contains(rc))
    setGeometry(rc);

  SetupScriptPath();

  connect(ui->pushButtonNext, SIGNAL(clicked(bool)), SLOT(OnButtonNext()));
  connect(ui->pushButtonPrev, SIGNAL(clicked(bool)), SLOT(OnButtonPrev()));
  connect(ui->pushButtonRegister, SIGNAL(clicked(bool)), SLOT(OnButtonRegister()));
  connect(ui->pushButtonClear, SIGNAL(clicked(bool)), SLOT(OnButtonClear()));

  setWindowTitle(WND_TITLE);

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
  s.setValue("MainWindow/Geometry", geometry());

  if (m_proc && m_proc->state() == QProcess::Running)
  {
    m_proc->kill();
    m_proc->waitForFinished();
  }
  delete ui;
}

void MainWindow::ShowDialog()
{
  hide();
  DialogWelcome dlg;
  dlg.setWindowTitle(WND_TITLE);
  if (dlg.exec() != QDialog::Accepted)
  {
    QTimer::singleShot(0, this, SLOT(close()));
    return;
  }

  DialogSelectFolder dlgSelect;
  dlgSelect.setWindowTitle(WND_TITLE);
  if (dlgSelect.exec() != QDialog::Accepted)
  {
    QTimer::singleShot(0, this, SLOT(close()));
    return;
  }

  show();
  m_nNumberOfExpectedPoints = dlgSelect.IsFourPoint()?4:2;
  ui->labelHeight->setVisible(m_nNumberOfExpectedPoints == 4);
  ui->lineEditRulerHeight->setVisible(m_nNumberOfExpectedPoints == 4);
  ui->widgetImageView->SetNumberOfExpectedPoints(m_nNumberOfExpectedPoints);
  ui->labelWidth->setText(m_nNumberOfExpectedPoints == 2? "Ruler Length (mm):":"Ruler Width (mm):");

  QString input_path = dlgSelect.GetInputPath();
  m_listInputFiles = QDir(input_path).entryInfoList(QDir::Files, QDir::Name);
  m_strOutputFolder = dlgSelect.GetOutputPath();

  UpdateIndex();

  if (!m_listInputFiles.isEmpty())
    LoadImage(m_nIndex);
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

  if (!QFile::exists(m_strPyScriptPath))
  {
    QMessageBox::critical(this, "Error", tr("Could not locate %1 script").arg(SCRIPT_FILENAME));
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }
}

void MainWindow::OnButtonNext()
{
  if (m_nIndex < m_listInputFiles.size()-1)
  {
    LoadImage(++m_nIndex);
    UpdateIndex();
  }
  else
  {
    if (QMessageBox::question(this, "", "No more images to register. Quit the program?") == QMessageBox::Yes)
      qApp->quit();
  }
}

void MainWindow::OnButtonPrev()
{
  LoadImage(--m_nIndex);
  UpdateIndex();
}

void MainWindow::OnButtonRegister()
{
  bool bOK;
  double dWidth = ui->lineEditRulerWidth->text().trimmed().toDouble(&bOK);
  double dHeight = 1;
  if (m_nNumberOfExpectedPoints == 4)
  {
    dHeight = ui->lineEditRulerHeight->text().trimmed().toDouble(&bOK);
    if (!bOK)
      dHeight = -1;
  }
  if (ui->widgetImageView->GetEditedPoints().size() < m_nNumberOfExpectedPoints)
    QMessageBox::warning(this, "", "Please click on the image to add another point");
  else if (!bOK || dWidth <= 0 || dHeight <= 0)
    QMessageBox::warning(this, "", "Please enter a valid value for ruler length");
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
    cmd << m_strPythonCmd << m_strPyScriptPath
        << "--in_img" << ui->widgetImageView->GetFilename()
        << "--points" << strList.join(" ")
        << "--width" << QString::number(dWidth)
        << "--height" << QString::number(dHeight)
        << "--out_dir" << m_strOutputFolder;
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
  ui->pushButtonRegister->setEnabled(false);
  ui->widgetImageView->ShowMessage("Processing...");
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->HideMessage();
  if (m_proc->exitStatus() == QProcess::NormalExit)
    ui->pushButtonNext->setEnabled(true);
  ui->pushButtonRegister->setEnabled(true);
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
}

void MainWindow::LoadImage(int n)
{
  QList<QPoint> pts;
  if (m_listData.size() > n)
    pts = m_listData[n];
  ui->widgetImageView->LoadImage(m_listInputFiles[n].absoluteFilePath(), "", pts);
}

void MainWindow::UpdateIndex()
{
  ui->labelIndex->setText(tr("%1 / %2").arg(m_nIndex+1).arg(m_listInputFiles.size()));
  ui->pushButtonPrev->setEnabled(m_nIndex > 0);
  ui->pushButtonNext->setEnabled(m_nIndex < m_listData.size());
}
