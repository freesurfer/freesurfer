/**
 * @file  DialogNewPointSet.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:46 $
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
#ifndef DIALOGNEWPOINTSET_H
#define DIALOGNEWPOINTSET_H

#include <QDialog>

namespace Ui
{
class DialogNewPointSet;
}

class LayerMRI;

class DialogNewPointSet : public QDialog
{
  Q_OBJECT

public:
  explicit DialogNewPointSet(QWidget *parent = 0);
  ~DialogNewPointSet();

  QString GetPointSetName();
  void SetPointSetName( const QString& name );

  LayerMRI* GetTemplate();

  int GetType();

protected slots:
  void OnOK();

private:
  Ui::DialogNewPointSet *ui;
};

#endif // DIALOGNEWPOINTSET_H
