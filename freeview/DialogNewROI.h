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
#ifndef DIALOGNEWROI_H
#define DIALOGNEWROI_H

#include <QDialog>

namespace Ui
{
class DialogNewROI;
}

class LayerMRI;

class DialogNewROI : public QDialog
{
  Q_OBJECT

public:
  explicit DialogNewROI(QWidget *parent = 0);
  ~DialogNewROI();

  QString GetROIName();
  void SetROIName( const QString& name );

  LayerMRI* GetTemplate();

protected slots:
  void OnOK();

private:
  Ui::DialogNewROI *ui;
};

#endif // DIALOGNEWROI_H
