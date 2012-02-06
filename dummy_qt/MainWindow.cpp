#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QProcess>
#include <QScrollBar>
#include <QDateTime>
#include <QFileInfo>
#include <QSettings>
#include <QCloseEvent>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_process = new QProcess(this);
    connect(m_process, SIGNAL(readyReadStandardOutput()),
            this, SLOT(OnStandardOutput()));
    connect(m_process, SIGNAL(readyReadStandardError()),
            this, SLOT(OnStandardOutput()));
    connect(m_process, SIGNAL(stateChanged(QProcess::ProcessState)),
            this, SLOT(UpdateStatus()));

    QSettings s;
    m_strLastDir = s.value("Application/LastDirectory").toString();
    restoreGeometry(s.value("Application/WindowGeometry").toByteArray());
}

MainWindow::~MainWindow()
{
    QSettings s;
    s.setValue("Application/LastDirectory", m_strLastDir);
    s.setValue("Application/WindowGeometry", saveGeometry());
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent * e)
{
  if (m_process->state() != QProcess::NotRunning && OnButtonAbort() == false)
  {
    e->ignore();
  }
  else
  {
    m_process->kill();
    e->accept();
  }
}

void MainWindow::OnButtonOpenInput()
{
  QString fn = QFileDialog::getOpenFileName(this, "Open", m_strLastDir);
  if (!fn.isEmpty())
  {
    ui->lineEditInput->setText(fn);
    ui->lineEditInput->setCursorPosition( ui->lineEditInput->text().size() );
    m_strLastDir = QFileInfo(fn).absolutePath();
  }
}

void MainWindow::OnButtonOpenOutput()
{
  QString fn = QFileDialog::getSaveFileName(this, "Save", m_strLastDir);
  if (!fn.isEmpty())
  {
    ui->lineEditOutput->setText(fn);
    ui->lineEditOutput->setCursorPosition( ui->lineEditOutput->text().size() );
  }
}

// Run the actual command here in a separate process
void MainWindow::OnButtonRun()
{
  QString cmd = QString("mri_convert %1 %2").arg(ui->lineEditInput->text().trimmed())
                .arg(ui->lineEditOutput->text().trimmed());
  m_process->start(cmd);

  AddMessage(QString("\n-----------%1-----------\n").arg(QDateTime::currentDateTime().toString(Qt::DefaultLocaleShortDate)),
             Highlight);
}

bool MainWindow::OnButtonAbort()
{
  if (QMessageBox::question(this, "Abort", "Do you really want to abort current process?",
                        QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
  {
    m_process->kill();
    AddMessage("Aborted by user", Error);
    return true;
  }
  else
    return false;
}

void MainWindow::UpdateStatus()
{
  bool ready = !ui->lineEditInput->text().trimmed().isEmpty() &&
      !ui->lineEditOutput->text().trimmed().isEmpty();
  bool busy = (m_process->state() != QProcess::NotRunning);
  ui->pushButtonRun->setEnabled(ready && !busy);
  ui->pushButtonAbort->setEnabled(ready && busy);
  this->setWindowTitle(busy?"mri_convert [running...]":"mri_convert");
}

void MainWindow::OnStandardOutput()
{
  AddMessage( m_process->readAllStandardOutput() );
  AddMessage( m_process->readAllStandardError(), Error);
}

void MainWindow::AddMessage(const QString& msg, int type)
{
  if (msg.trimmed().isEmpty())
    return;

  QString strg = msg;
  QString color = "black";
  switch (type)
  {
  case Highlight:
    color = "green";
    break;
  case Error:
    color = "#C00";
    break;
  case Warning:
    color = "#CC0";
    break;
  }
  if (type != Normal)
    strg = QString("<span style=\"font-weight:bold; color:%1;\">%2</span>").arg(color).arg(msg);
  else
    strg = QString("<span style=\"font-weight:normal; color:%1;\">%2</span>").arg(color).arg(msg);
  strg.replace("\n", "<br />");
  ui->textEditLog->append(strg);
  QScrollBar* vbar = ui->textEditLog->verticalScrollBar();
  if (!vbar->isSliderDown())
    vbar->setValue(vbar->maximum());
}
