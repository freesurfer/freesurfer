/**
 * @file  DialogPreferences.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.16 $
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
#ifndef DIALOGPREFERENCES_H
#define DIALOGPREFERENCES_H

#include <QDialog>
#include <QVariantMap>
#include "UIUpdateHelper.h"

namespace Ui
{
class DialogPreferences;
}

class QAbstractButton;

class DialogPreferences : public QDialog, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit DialogPreferences(QWidget *parent = 0);
  ~DialogPreferences();

  void SetSettings(const QVariantMap& map);

  QVariantMap GetSettings();

protected slots:
  void OnClicked(QAbstractButton* btn);

private:
  Ui::DialogPreferences *ui;
};

#endif // DIALOGPREFERENCES_H
