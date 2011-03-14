/**
 * @file  DialogSavePointSet.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

  void SetFileName(const QString& fn);

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

private:
  Ui::DialogSavePointSet *ui;
  QString m_strLastDir;
};

#endif // DIALOGSAVEPOINTSET_H
