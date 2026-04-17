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
#include "CommonDef.h"
#include <QTemporaryDir>
#include <QScreen>
#include <QPainter>

#define WND_TITLE "Mask Extraction"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent), m_bProfiling(false), ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);

  m_strPythonCmd = PY_COMMAND;

  // setup script
  static QTemporaryDir dir;
  m_strTempFolder = dir.path();
//  m_strTempFolder = "/Users/rpwang/test_data/photo_recon/tmp";
  m_strPyScriptPath = QFileInfo(dir.path(),"func_mask_extraction.py").absoluteFilePath();
  QFile::copy(":/func_mask_extraction.py", m_strPyScriptPath);

  if (!QFile::exists(m_strPyScriptPath))
  {
    QMessageBox::critical(this, "Error", "Could not locate func_mask_extraction.py script");
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }

  m_listStockColors << QColor(255,100,100) << QColor(255,255,100) << QColor(100,255,100)
                    << QColor(110,245,255) << QColor(75,100,255) << QColor(255,128,0)
                    << QColor(100,150,170) << QColor(120,60,220);

  ui->widgetImageView->SetEditMode(WidgetImageView::EM_SELECT_MASK);
  ui->sliderMaskOpacity->setValue(ui->widgetImageView->GetMaskOpacity()*100);

  connect(ui->pushButtonNext, SIGNAL(clicked(bool)), SLOT(OnButtonNext()));
  connect(ui->pushButtonPrev, SIGNAL(clicked(bool)), SLOT(OnButtonPrev()));
  connect(ui->pushButtonCreateMask, SIGNAL(clicked(bool)), SLOT(OnButtonCreateMask()));
  connect(ui->pushButtonClear, SIGNAL(clicked(bool)), SLOT(OnButtonClear()));
  connect(ui->sliderMaskOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderMaskOpacity(int)));

  setWindowTitle(WND_TITLE);
  QSettings s;
  QRect rc = s.value("MainWindow/Geometry").toRect();
  if (rc.isValid() && QGuiApplication::primaryScreen()->geometry().contains(rc))
    setGeometry(rc);

  m_proc = new QProcess(this);
  connect(m_proc, SIGNAL(started()), SLOT(OnProcessStarted()));
  connect(m_proc, SIGNAL(readyReadStandardOutput()), SLOT(OnProcessOutputMessage()));
  connect(m_proc, SIGNAL(readyReadStandardError()), SLOT(OnProcessErrorMessage()));
  connect(m_proc, SIGNAL(finished(int, QProcess::ExitStatus)), SLOT(OnProcessFinished()));
  connect(m_proc, SIGNAL(errorOccurred(QProcess::ProcessError)), SLOT(OnProcessError(QProcess::ProcessError)));
}

MainWindow::~MainWindow()
{
  QSettings s;
  s.setValue("MainWindow/Geometry", geometry());

  if (m_proc && m_proc->state() == QProcess::Running)
  {
//    m_proc->kill();
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

  QString path = dlgSelect.GetInputPath();
  m_listInputFiles = QDir(path).entryInfoList(QDir::Files, QDir::Name);

  m_strOutputFolder = dlgSelect.GetOutputPath();

  QStringList cmd;
  cmd << m_strPythonCmd << m_strPyScriptPath
       << "--input" << QString("\"%1\"").arg(m_listInputFiles.first().absoluteFilePath())
      << "--output" << QString("\"%1\"").arg(m_strTempFolder)
      << "--checkpoint /Users/rpwang/test_data/segment-anything/sam_vit_h_4b8939.pth --model-type default --device cpu";
  m_proc->start(cmd.join(" "));
//  qDebug() << cmd.join(" ");

  UpdateIndex();
}

void MainWindow::OnButtonNext()
{
  if (m_nIndex < m_listInputFiles.size()-1)
  {
    ui->pushButtonNext->setEnabled(false);
    LoadImage(++m_nIndex);
    UpdateIndex();
  }
  else
  {
    if (QMessageBox::question(this, "", "No more images to process. Quit the program?") == QMessageBox::Yes)
      qApp->quit();
  }
}

void MainWindow::OnButtonPrev()
{
  LoadImage(--m_nIndex);
  UpdateIndex();
}

void MainWindow::OnButtonCreateMask()
{
  QList<QImage> list = ui->widgetImageView->GetSelectedMasks();
  if (!list.isEmpty())
  {
    QImage img(list[0].size(), QImage::Format_ARGB32);
    img.fill(QColor(0,0,0,0));
    QPainter p(&img);
    foreach (QImage im, list)
      p.drawImage(0, 0, im);
    p.end();
    img.save(m_strOutputFolder + "/" + m_listInputFiles[m_nIndex].completeBaseName()+"_mask.png");
    ui->pushButtonNext->setEnabled(true);
  }
}

void MainWindow::OnButtonClear()
{
  ui->widgetImageView->ClearEdits();
  ui->pushButtonNext->setEnabled(false);
}

void MainWindow::LoadImage(int n)
{
  QString fn = m_listInputFiles[n].absoluteFilePath();
  QFileInfoList info_list = QDir(m_strTempFolder + "/" + m_listInputFiles[n].completeBaseName()).entryInfoList(QDir::Files, QDir::Name);
  QStringList list;
  foreach (QFileInfo fi, info_list)
    list << fi.absoluteFilePath();
//  qDebug() << fn << list;
  ui->widgetImageView->LoadImage(fn, list);
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
  ui->widgetImageView->ShowMessage("Segmenting...");
  ui->widgetImageView->setDisabled(true);
  ui->pushButtonCreateMask->setEnabled(false);
  if (m_bProfiling)
    m_timer.start();
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->setEnabled(true);
  if (!m_listInputFiles.isEmpty())
    LoadImage(0);

  ui->widgetImageView->HideMessage();
  ui->pushButtonCreateMask->setEnabled(true);
  if (m_bProfiling)
  {
    qDebug() << "Elapsed time by image processing: " << m_timer.elapsed() << "ms";
    m_timer.invalidate();
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

void MainWindow::UpdateIndex()
{
  ui->labelIndex->setText(tr("%1 / %2").arg(m_nIndex+1).arg(m_listInputFiles.size()));
  ui->pushButtonPrev->setEnabled(m_nIndex > 0 && m_proc->state() != QProcess::Running);
  ui->pushButtonNext->setEnabled(m_nIndex < m_listInputFiles.size());
}

void MainWindow::OnSliderMaskOpacity(int val)
{
  ui->widgetImageView->SetMaskOpacity(1.0*val/100);
}
