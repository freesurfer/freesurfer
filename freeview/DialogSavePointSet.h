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
#ifndef DIALOGSAVEPOINTSET_H
#define DIALOGSAVEPOINTSET_H

#include <QDialog>

namespace Ui
{
class DialogSavePointSet;
}

class DialogSavePointSet : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSavePointSet(QWidget *parent = 0);
  ~DialogSavePointSet();

  void SetFileName(const QString& fn, int type);

  void SetLastDir(const QString& dir)
  {
    m_strLastDir = dir;
  }

  void SetType( int nType );

  int GetType();

  QString GetFileName();

protected slots:
  void OnOK();
  void OnOpen();
  void OnRadioButtonFileType();

private:
  Ui::DialogSavePointSet *ui;
  QString m_strLastDir;

  bool  m_bRemind;
};

#endif // DIALOGSAVEPOINTSET_H
