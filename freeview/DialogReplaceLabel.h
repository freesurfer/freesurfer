/**
 * @brief Dialog window to replace label.
 *
 */
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
 *
 */

#ifndef DIALOGREPLACELABEL_H
#define DIALOGREPLACELABEL_H

#include <QDialog>

namespace Ui {
  class DialogReplaceLabel;
}

class LayerMRI;

class DialogReplaceLabel : public QDialog
{
  Q_OBJECT

public:
  explicit DialogReplaceLabel(QWidget *parent = 0);
  ~DialogReplaceLabel();

  double GetOriginalValue();

  double GetNewValue();

  bool ReplaceSingleSlice();

public slots:
  void OnReplace();

private:
  Ui::DialogReplaceLabel *ui;
};

#endif // DIALOGREPLACELABEL_H
