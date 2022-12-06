#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "CommonDef.h"
#include <QTemporaryDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QSettings>
#include <QDebug>
#include <QDesktopWidget>
#include <QFileDialog>
#include "DialogWelcome.h"
#include <QTimer>

#define WND_TITLE "Calibration With Fiducials"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  m_strPythonCmd = PY_COMMAND;

  connect(ui->pushButtonClear, SIGNAL(clicked(bool)), SLOT(OnButtonClear()));
  connect(ui->pushButtonGo, SIGNAL(clicked(bool)), SLOT(OnButtonGo()));
  connect(ui->pushButtonLoadImage, SIGNAL(clicked(bool)), SLOT(OnButtonLoad()));

  // setup script
  static QTemporaryDir dir;
  m_strTempFolder = dir.path();
  m_strPyScriptPath = QFileInfo(dir.path(),"func_fiducials_calibration.py").absoluteFilePath();
  QFile::copy(":/func_fiducials_calibration.py", m_strPyScriptPath);

  if (!QFile::exists(m_strPyScriptPath))
  {
    QMessageBox::critical(this, "Error", "Could not locate func_fiducials_calibration.py script");
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }

  ui->widgetImageView->SetEditMode(WidgetImageView::EM_CALIBRATION);
  connect(ui->widgetImageView, SIGNAL(CalibrationReady(QList<QPoint>)), SLOT(OnCalibrationReady(QList<QPoint>)));

  setWindowTitle(WND_TITLE);

  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  m_strLastDir = s.value("LastDir").toString();
  if (rc.isValid() && QApplication::desktop()->screenGeometry(this).contains(rc))
    setGeometry(rc);

  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()));
  connect(m_proc, SIGNAL(error(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("MainWindow/Geometry", geometry());
  s.setValue("LastDir", m_strLastDir);

  if (m_proc && m_proc->state() == QProcess::Running)
  {
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

  show();
  OnButtonLoad();
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
  ui->widgetImageView->ShowMessage(str);
}

void MainWindow::OnProcessStarted()
{
  ui->widgetImageView->setDisabled(true);
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->setEnabled(true);
  ui->widgetImageView->HideMessage();

  QMessageBox dlg(QMessageBox::Information, "Calibration", "Calibration finished. Click Exit to quit.", QMessageBox::Ok, this);
  dlg.setButtonText(QMessageBox::Ok, "Exit");
  dlg.exec();
  QTimer::singleShot(0, qApp, SLOT(quit()));
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
  ui->pushButtonGo->setEnabled(false);
}

void MainWindow::OnButtonLoad()
{
  QString fn = QFileDialog::getOpenFileName(this, "Select Image", m_strLastDir);
  if (!fn.isEmpty())
  {
    ui->widgetImageView->LoadImage(fn);
    ui->widgetImageView->ClearEdits();
    ui->pushButtonGo->setEnabled(false);
    m_strLastDir = QFileInfo(fn).absolutePath();
  }
}

void MainWindow::OnButtonGo()
{
  double w = 1, h = 1;
  bool bOK;
  w = ui->lineEditWidth->text().trimmed().toDouble(&bOK);
  if (bOK)
    h = ui->lineEditHeight->text().trimmed().toDouble(&bOK);
  if (!bOK)
  {
    QMessageBox::warning(this, "Error", "Please enter valid width and height");
    return;
  }

  QString fn = QFileDialog::getSaveFileName(this, "Specify Output npz File", m_strLastDir, "*.npz");
  if (!fn.isEmpty())
  {
    QStringList cmd;
    QStringList list;
    foreach (QPoint pt, m_listPoints)
      list << QString::number(pt.x()) << QString::number(pt.y());
    cmd << m_strPythonCmd << m_strPyScriptPath
        << "--in_img" << ui->widgetImageView->GetFilename()
        << "--points" << list.join(" ")
        << "--width" << QString::number(w)
        << "--height" << QString::number(h)
        << "--out_file" << fn;
    m_proc->start(cmd.join(" "));
    ui->widgetImageView->ShowMessage("Processing...");
  }
}

void MainWindow::OnCalibrationReady(const QList<QPoint>& pts)
{
  m_listPoints = pts;
  ui->pushButtonGo->setEnabled(true);
}
