/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "TermWidget.h"
#include "ui_TermWidget.h"
#include "MainWindow.h"
#include "MyCmdLineParser.h"
#include <QTimer>
#include <QDebug>
#include <QSettings>
#include <QScrollBar>
#include <iostream>
#include <stdlib.h>

#define BUFFER_SIZE 8192
char stdout_buffer[BUFFER_SIZE] = {0};
char stderr_buffer[BUFFER_SIZE] = {0};

void outcallback( const char* ptr, std::streamsize count, void* pBuffer )
{
  Q_UNUSED(count)
  QString* str = static_cast<QString*>(pBuffer);
  *str += ptr;
}

void errcallback( const char* ptr, std::streamsize count, void* pBuffer )
{
  Q_UNUSED(count)
  QString* str = static_cast<QString*>(pBuffer);
  *str += ptr;
}

QStringList Known_Shell_Cmds;

TermWidget::TermWidget(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::TermWidget),
  m_stdOut(0),
  m_stdErr(0),
  m_stdinNotifier(NULL)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Tool );
  SetDarkTheme(true);
  QFont font = ui->textCommand->font();
  font.setFamily("Menlo");
  font.setPixelSize(12);
  font.setStyleHint(QFont::TypeWriter);
  ui->textCommand->setFont(font);
  ui->textLog->setFont(font);
  ui->textLog->setWordWrapMode(QTextOption::WrapAnywhere);
  ui->textCommand->setFocus();
  connect(ui->textCommand, SIGNAL(CommandTriggered(QString)),
          this, SLOT(OnCommandTriggered(QString)), Qt::QueuedConnection);

  QSettings settings;
  this->restoreGeometry(settings.value("/CommandConsole/Geometry").toByteArray());

  SetRedirectStdOutput(true);
  SetDuplicateStdOutput(true);

  Known_Shell_Cmds << "ls" << "pwd" << "cd" << "cp" << "dir" << "copy";

  QTimer* timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(OnTimeOut()));
  timer->start(200);
}

TermWidget::~TermWidget()
{
  QSettings settings;
  settings.setValue("/CommandConsole/Geometry", this->saveGeometry());

  SetRedirectStdOutput(false);
  delete ui;
}

void TermWidget::EnableListeningStdin()
{
  if (!m_stdinNotifier)
  {
    m_stdinNotifier = new QSocketNotifier(STDIN_FILENO, QSocketNotifier::Read, this);
    connect(m_stdinNotifier, SIGNAL(activated(int)), this, SLOT(OnStdinActivated()));
  }
}

void TermWidget::SetRedirectStdOutput(bool bRedir)
{
  m_bRedirectStdOutput = bRedir;
  if (bRedir)
  {
    if (m_stdOut)
    {
      delete m_stdOut;
    }
    if (m_stdErr)
    {
      delete m_stdErr;
    }

    m_stdOut = new StdRedirector<>( std::cout, outcallback, &m_bufferStdOut );
    m_stdErr = new StdRedirector<>( std::cerr, errcallback, &m_bufferStdErr );
    fflush(stdout);
    setvbuf(stdout, stdout_buffer, _IOLBF, BUFFER_SIZE);
    fflush(stderr);
    setvbuf(stderr, stderr_buffer, _IOLBF, BUFFER_SIZE);
    stdout_buffer[0] = 0;
    stderr_buffer[0] = 0;
  }
  else
  {
    delete m_stdOut;
    delete m_stdErr;
    m_stdOut = 0;
    m_stdErr = 0;
    fflush(stdout);
    fflush(stderr);
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
  }
}

void TermWidget::OnCommandTriggered(const QString &cmd)
{
  QString strg = cmd;
  if (strg.trimmed().isEmpty())
  {
    return;
  }

  if (strg == "abort")
  {
    MainWindow::GetMainWindow()->AbortScripts();
    AppendErrorString("Script aborted");
    return;
  }

  if (strg[strg.size()-1] == '&')
  {
    strg.resize(strg.size()-1);
  }
  ui->textLog->append(QString("<span style=\"color:%1;\">%2</span>")
                      .arg(m_strLogColor)
                      .arg(ui->textCommand->GetPrompt() + strg));
  ScrollToBottom();
  m_bufferStdOut.clear();
  m_bufferStdErr.clear();
  QStringList args = cmd.split(" ", QString::SkipEmptyParts);
  if (Known_Shell_Cmds.contains(args[0].toLower()))
  {
    AppendErrorString("This is not a shell. Only freeview commands are supported. Type '-h' for all the available commands.\n");
  }
  else if (args[0].toLower() == "clear")
  {
    ui->textLog->clear();
  }
  else
  {
    if (strg[0] == '-')
      MainWindow::GetMainWindow()->ParseCommand(QString("freeview ") + strg);
    else
      MainWindow::GetMainWindow()->AddScript(strg.split(" ", QString::SkipEmptyParts));

    if ( MainWindow::GetMainWindow()->IsBusy())
    {
      AppendErrorString("Still busy. Command is added to queue and will be executed later.\n");
    }
  }
}

void TermWidget::OnStdinActivated()
{
  if (!m_stdinNotifier)
    return;

  m_stdinNotifier->setEnabled(false);

  char sbuf[8192];
  fgets(sbuf, sizeof(sbuf), stdin);
  fflush(stdin);

  QString line = QString::fromUtf8(sbuf).trimmed();
  if (line.indexOf("freeview") == 0)
  {
    line = line.mid(8).trimmed();
    OnCommandTriggered(line);
  }
  m_stdinNotifier->setEnabled(true);
}

void TermWidget::OnTimeOut()
{
  if (!m_bRedirectStdOutput)
  {
    return;
  }

  if (stdout_buffer[0] != 0)
  {
    m_bufferStdOut += stdout_buffer;
    fflush(stdout);
    stdout_buffer[0] = 0;
  }
  if (!m_bufferStdOut.isEmpty())
  {
    ui->textLog->append(m_bufferStdOut);
    ScrollToBottom();
    if (m_bDuplicateStdOutput)
    {
      delete m_stdOut;
      setbuf(stdout, NULL);
      std::cout << qPrintable(m_bufferStdOut);
      fflush(0);
      m_bufferStdOut.clear();
      m_stdOut = new StdRedirector<>( std::cout, outcallback, &m_bufferStdOut );
      setvbuf(stdout, stdout_buffer, _IOLBF, BUFFER_SIZE);
    }
    else
    {
      m_bufferStdOut.clear();
    }
  }

  if (stderr_buffer[0] != 0)
  {
    m_bufferStdErr += stderr_buffer;
    fflush(stderr);
    stderr_buffer[0] = 0;
  }
  if (!m_bufferStdErr.isEmpty())
  {
    AppendErrorString(m_bufferStdErr);
    ScrollToBottom();
    if (m_bDuplicateStdOutput)
    {
      delete m_stdErr;
      setbuf(stderr, NULL);
      std::cerr << qPrintable(m_bufferStdErr);
      fflush(0);
      m_bufferStdErr.clear();
      m_stdErr = new StdRedirector<>( std::cerr, outcallback, &m_bufferStdErr );
      setvbuf(stderr, stderr_buffer, _IOLBF, BUFFER_SIZE);
    }
    else
    {
      m_bufferStdErr.clear();
    }
    MainWindow::GetMainWindow()->SetHadError(true);
  }
}

void TermWidget::AppendErrorString(const QString &strg_in)
{
  QString strg = strg_in;
  strg.replace("\n", "<br />");
  ui->textLog->append(QString("<span style=\"color:%1;\">%2</span>")
                      .arg(m_strErrorColor)
                      .arg(strg.trimmed()));
}

void TermWidget::SetDarkTheme(bool bDark)
{
  QScrollBar* sb = ui->textLog->verticalScrollBar();
  if (bDark)
  {
    QPalette pal = ui->textLog->palette();
    pal.setColor(QPalette::Text, QColor(0x909090));
    pal.setColor(QPalette::Base, QColor(0x333333));
    ui->textLog->setPalette(pal);
    m_strLogColor = "#c0c0c0";
    m_strErrorColor = "#ee4444";
    QString text = ui->textLog->toHtml();
    text.replace("#444444", m_strLogColor);
    text.replace("#ff0000", m_strErrorColor);
    int pos = sb->sliderPosition();
    ui->textLog->setHtml(text);
    sb->setSliderPosition(pos);

    pal = ui->textCommand->palette();
    pal.setColor(QPalette::Text, QColor(0xd0d0d0));
    pal.setColor(QPalette::Base, QColor(0x222222));
    ui->textCommand->setPalette(pal);
  }
  else
  {
    QPalette pal = ui->textLog->palette();
    pal.setColor(QPalette::Text, QColor(0x777777));
    pal.setColor(QPalette::Base, QColor(0xeeeeee));
    ui->textLog->setPalette(pal);
    m_strLogColor = "#444444";
    m_strErrorColor = "#ff0000";
    QString text = ui->textLog->toHtml();
    text.replace("#c0c0c0", m_strLogColor);
    text.replace("#ee4444", m_strErrorColor);
    int pos = sb->sliderPosition();
    ui->textLog->setHtml(text);
    sb->setSliderPosition(pos);

    pal = ui->textCommand->palette();
    pal.setColor(QPalette::Text, QColor(0x010101));
    pal.setColor(QPalette::Base, QColor(0xffffff));
    ui->textCommand->setPalette(pal);
  }
}

void TermWidget::ScrollToBottom()
{
  QScrollBar* sb = ui->textLog->verticalScrollBar();
  if (sb)
  {
    sb->setSliderPosition(sb->maximum());
  }
}

bool TermWidget::GetDarkTheme()
{
  return (m_strLogColor != "#444444");
}
