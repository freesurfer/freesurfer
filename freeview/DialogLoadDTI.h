/**
 * @file  DialogLoadDTI.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/09/25 18:45:12 $
 *    $Revision: 1.15 $
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
#ifndef DIALOGLOADDTI_H
#define DIALOGLOADDTI_H

#include <QDialog>

namespace Ui
{
class DialogLoadDTI;
}

class DialogLoadDTI : public QDialog
{
  Q_OBJECT

public:
  explicit DialogLoadDTI(QWidget *parent = 0);
  ~DialogLoadDTI();

  QString GetVectorFileName();
  QString GetFAFileName();
  QString GetRegFileName();
  QString GetEigenvalueFileName();

  bool IsToResample();

  void SetLastDir( const QString& dir )
  {
    m_strLastDir = dir;
  }

  void Initialize( bool bResample, bool bEnableCheckBox );

protected slots:
  void OnOK();
  void OnButtonVector();
  void OnButtonFA();
  void OnButtonRegistration();
  void OnButtonEigenvalue();

private:
  Ui::DialogLoadDTI *ui;
  QString m_strLastDir;
};

#endif // DIALOGLOADDTI_H
