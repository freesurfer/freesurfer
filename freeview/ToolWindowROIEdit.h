/**
 * @file  ToolWindowROIEdit.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:18 $
 *    $Revision: 1.2 $
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

#ifndef TOOLWINDOWROIEDIT_H
#define TOOLWINDOWROIEDIT_H

#include <QWidget>

namespace Ui
{
class ToolWindowROIEdit;
}

class ToolWindowROIEdit : public QWidget
{
  Q_OBJECT

public:
  explicit ToolWindowROIEdit(QWidget *parent = 0);
  ~ToolWindowROIEdit();

  void UpdateWidgets();

protected:
  virtual void showEvent(QShowEvent *);

protected slots:
  void OnEditMode(QAction *act);

private:
  Ui::ToolWindowROIEdit *ui;
};

#endif // TOOLWINDOWROIEDIT_H
