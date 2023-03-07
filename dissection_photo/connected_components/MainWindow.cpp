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

#define WND_TITLE "Connected Components"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent), m_bProfiling(false), ui(new Ui::MainWindow), m_nIndex(0)
{
  ui->setupUi(this);

  m_strPythonCmd = PY_COMMAND;

  // setup script
  static QTemporaryDir dir;
  m_strTempFolder = dir.path();
  m_strPyScriptPath = QFileInfo(dir.path(),"func_mask_to_cc.py").absoluteFilePath();
  QFile::copy(":/func_mask_to_cc.py", m_strPyScriptPath);

  if (!QFile::exists(m_strPyScriptPath))
  {
    QMessageBox::critical(this, "Error", "Could not locate func_mask_to_cc.py script");
    QTimer::singleShot(0, qApp, SLOT(quit()));
  }

  m_listStockColors << QColor(255,100,100) << QColor(255,255,100) << QColor(100,255,100)
                    << QColor(110,245,255) << QColor(75,100,255) << QColor(255,128,0)
                    << QColor(100,150,170) << QColor(120,60,220);

  ui->widgetImageView->SetEditMode(WidgetImageView::EM_REGION);
  connect(ui->widgetImageView, SIGNAL(LastRegionEdited(int)),
                                      SLOT(OnLastRegionEdited(int)), Qt::QueuedConnection);

  connect(ui->pushButtonNext, SIGNAL(clicked(bool)), SLOT(OnButtonNext()));
  connect(ui->pushButtonPrev, SIGNAL(clicked(bool)), SLOT(OnButtonPrev()));
  connect(ui->pushButtonCreateMask, SIGNAL(clicked(bool)), SLOT(OnButtonCreateMask()));
  connect(ui->pushButtonClear, SIGNAL(clicked(bool)), SLOT(OnButtonClear()));

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

  path = dlgSelect.GetMaskPath();
  m_listMaskFiles = QDir(path).entryInfoList(QDir::Files, QDir::Name);

  m_strOutputFolder = dlgSelect.GetOutputPath();

  QStringList cmd;
  cmd << m_strPythonCmd << m_strPyScriptPath
       << "--in_dir" << QString("\"%1\"").arg(dlgSelect.GetMaskPath())
       << "--out_dir" << QString("\"%1\"").arg(m_strTempFolder);
  m_proc->start(cmd.join(" "));

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
  QList<RECT_REGION> list = ui->widgetImageView->GetEditedRegions();
  if (m_nIndex < m_listData.size())
    m_listData[m_nIndex] = list;
  else
    m_listData << list;

  QString fn = QFileInfo(ui->widgetImageView->GetFilename()).fileName();
  fn.replace(QString(".")+QFileInfo(fn).suffix(), "_mask.npy", Qt::CaseInsensitive);
  m_maskProcessor.SaveToNpy(QFileInfo(m_strOutputFolder, fn).absoluteFilePath());
  ui->pushButtonNext->setEnabled(true);
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

void MainWindow::OnButtonClear()
{
  ui->widgetImageView->ClearEdits();
  ui->pushButtonNext->setEnabled(false);
  m_maskProcessor.ClearBuffer();
}

void MainWindow::LoadImage(int n)
{
  QList<RECT_REGION> rects;
  if (m_listData.size() > n)
    rects = m_listData[n];
  ui->widgetImageView->LoadImage(m_listInputFiles[n].absoluteFilePath(), m_listMaskFiles[n].absoluteFilePath(),
                                 QList<QPoint>(), rects);
  QString mask_fn = m_listMaskFiles[n].fileName();
  mask_fn.replace(QString(".")+QFileInfo(mask_fn).suffix(), ".npz");
  if (!m_maskProcessor.Load(QFileInfo(m_strTempFolder, mask_fn).absoluteFilePath()))
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
  ui->widgetImageView->ShowMessage("Initializing...");
  ui->widgetImageView->setDisabled(true);
  ui->pushButtonCreateMask->setEnabled(false);
  if (m_bProfiling)
    m_timer.start();
}

void MainWindow::OnProcessFinished()
{
  ui->widgetImageView->setEnabled(true);
  if (!m_listInputFiles.isEmpty() && !m_listMaskFiles.isEmpty())
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
  ui->pushButtonNext->setEnabled(m_nIndex < m_listData.size());
}
