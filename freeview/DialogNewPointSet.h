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
#ifndef DIALOGNEWPOINTSET_H
#define DIALOGNEWPOINTSET_H

#include <QDialog>

namespace Ui
{
class DialogNewPointSet;
}

class LayerMRI;
class LayerSurface;

class DialogNewPointSet : public QDialog
{
  Q_OBJECT

public:
  explicit DialogNewPointSet(QWidget *parent = 0);
  ~DialogNewPointSet();

  QString GetPointSetName();
  void SetPointSetName( const QString& name );

  LayerMRI* GetTemplate();
  LayerSurface* GetTemplateSurface();

  int GetType();

protected slots:
  void OnOK();

private:
  Ui::DialogNewPointSet *ui;
};

#endif // DIALOGNEWPOINTSET_H
