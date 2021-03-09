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
#ifndef COMMANDEDIT_H
#define COMMANDEDIT_H

#include <QPlainTextEdit>
#include <QStringList>

class CommandEdit : public QTextEdit
{
  Q_OBJECT
public:
  explicit CommandEdit(QWidget *parent = 0);
  ~CommandEdit();

  QString GetPrompt()
  {
    return m_strPrompt;
  }

signals:
  void CommandTriggered(const QString& cmd);

protected slots:
  void OnSelectionChanged();

protected:
  virtual void keyPressEvent(QKeyEvent *e);
  virtual void mousePressEvent(QMouseEvent *e);
  virtual void mouseMoveEvent(QMouseEvent *e);

private:
  void ProcessCommandInput();
  QStringList m_listHistory;
  int         m_nPosInHistory;
  QString     m_strPrompt;
  QString     m_strTempCommand;
};

#endif // COMMANDEDIT_H
