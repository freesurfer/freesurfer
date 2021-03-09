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
#include "CommandEdit.h"
#include <QSettings>

#ifdef Q_OS_MAC
#define CONTROL_KEY Qt::MetaModifier
#else
#define CONTROL_KEY Qt::ControlModifier
#endif

CommandEdit::CommandEdit(QWidget *parent) :
  QTextEdit(parent)
{
  m_strPrompt = "[FreeView] ";
  this->append(m_strPrompt);
  this->setCursorWidth(7);
  this->setWordWrapMode(QTextOption::WrapAnywhere);
  connect(this, SIGNAL(selectionChanged()), this, SLOT(OnSelectionChanged()));
  QSettings settings;
  m_listHistory = settings.value("/CommandConsole/History").toStringList();
  m_nPosInHistory = m_listHistory.size();
  this->setAcceptRichText(false);
}

CommandEdit::~CommandEdit()
{
  QSettings settings;
  while (m_listHistory.size() > 50)
  {
    m_listHistory.removeFirst();
  }
  settings.setValue("/CommandConsole/History", m_listHistory);
}

void CommandEdit::keyPressEvent(QKeyEvent *e)
{
  int keyCode = e->key();
  if (keyCode == Qt::Key_Left || keyCode == Qt::Key_Backspace)
  {
    QTextCursor cursor = this->textCursor();
    if (cursor.position() <= m_strPrompt.size())
    {
      e->ignore();
      return;
    }
  }
  else if (keyCode == Qt::Key_Up)
  {
    if (m_nPosInHistory == m_listHistory.size())
    {
      m_strTempCommand = this->toPlainText().mid(m_strPrompt.size());
    }
    if (m_nPosInHistory > 0)
    {
      m_nPosInHistory--;
      this->clear();
      this->append(m_strPrompt + m_listHistory[m_nPosInHistory]);
    }
  }
  else if (keyCode == Qt::Key_Down)
  {
    if ( m_nPosInHistory <= m_listHistory.size()-1)
    {
      m_nPosInHistory++;
      this->clear();
      if (m_nPosInHistory < m_listHistory.size())
      {
        this->append(m_strPrompt + m_listHistory[m_nPosInHistory]);
      }
      else
      {
        this->append(m_strPrompt + m_strTempCommand);
      }
    }
  }
  else if (keyCode == Qt::Key_Home ||
           (keyCode == Qt::Key_A && e->modifiers() & CONTROL_KEY))
  {
    QTextCursor cursor = this->textCursor();
    cursor.setPosition(m_strPrompt.size(), QTextCursor::MoveAnchor);
    this->setTextCursor(cursor);
    e->ignore();
    return;
  }
  else if (keyCode == Qt::Key_C && e->modifiers() & CONTROL_KEY)
  {
    this->clear();
    this->append(m_strPrompt);
    e->ignore();
    return;
  }
  else if (keyCode == Qt::Key_Return)
  {
    ProcessCommandInput();
    e->ignore();
    return;
  }
  else if (keyCode == Qt::Key_Tab)
  {
    e->ignore();
    return;
  }
  QTextEdit::keyPressEvent(e);
}

void CommandEdit::mousePressEvent(QMouseEvent *e)
{
  QTextCursor cursor = this->cursorForPosition(e->pos());
  if (cursor.position() < m_strPrompt.size())
  {
    e->ignore();
    return;
  }
  QTextEdit::mousePressEvent(e);
}

void CommandEdit::mouseMoveEvent(QMouseEvent *e)
{
  QTextCursor cursor = this->cursorForPosition(e->pos());
  if (cursor.position() < m_strPrompt.size())
  {
    e->ignore();
    return;
  }
  QTextEdit::mouseMoveEvent(e);
}

void CommandEdit::ProcessCommandInput()
{
  m_strTempCommand.clear();
  QString strg = this->toPlainText();
  QStringList list = strg.mid(m_strPrompt.size()).split(QRegExp("\\s+"), QString::SkipEmptyParts);
  if (!list.isEmpty() && (list[0].toLower() == "freeview" || list[0].toLower() == "fv") )
  {
    list.removeFirst();
  }
  strg = list.join(" ");
  this->clear();
  this->append(m_strPrompt);
  if ( m_listHistory.isEmpty() || m_listHistory.last() != strg)
  {
    m_listHistory << strg;
  }
  m_nPosInHistory = m_listHistory.size();
  emit CommandTriggered(strg);
}

void CommandEdit::OnSelectionChanged()
{
  QTextCursor cursor = this->textCursor();
  bool bToMove = false;
  int nStart = cursor.selectionStart();
  int nEnd = cursor.selectionEnd();
  if ( nStart < m_strPrompt.size())
  {
    bToMove = true;
    nStart = m_strPrompt.size();
  }
  if ( nEnd < m_strPrompt.size())
  {
    bToMove = true;
    nEnd = m_strPrompt.size();
  }
  if (bToMove)
  {
    cursor.setPosition(nStart, QTextCursor::MoveAnchor);
    cursor.setPosition(nEnd, QTextCursor::KeepAnchor);
    this->setTextCursor(cursor);
  }
}
