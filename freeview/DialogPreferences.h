/*
 * Original Author: Ruopeng Wang
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
class QComboBox;

class DialogPreferences : public QDialog, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit DialogPreferences(QWidget *parent = 0);
  ~DialogPreferences();

  void SetSettings(const QVariantMap& map);

  QVariantMap GetSettings();
  static void SetActionShortcut(QAction* act, const QString& text);

protected slots:
  void OnClicked(QAbstractButton* btn);
  void OnComboShortcutChanged(const QString& text);

private:
  void SetCurrentComboText(QComboBox* combo, const QString& text);

  Ui::DialogPreferences *ui;
};

#endif // DIALOGPREFERENCES_H
