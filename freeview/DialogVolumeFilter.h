/**
 * @file  DialogVolumeFilter.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.5 $
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

#ifndef DIALOGVOLUMEFILTER_H
#define DIALOGVOLUMEFILTER_H

#include <QDialog>

namespace Ui
{
class DialogVolumeFilter;
}

class VolumeFilter;

class DialogVolumeFilter : public QDialog
{
  Q_OBJECT

public:
  explicit DialogVolumeFilter(QWidget *parent = 0);
  ~DialogVolumeFilter();

  void SetFilter( VolumeFilter* filter );

  int GetKernelSize();

  double GetSigma();

  void SetSigma( double dvalue );

  void ShowSigma( bool bShow );

protected slots:
  void OnOK();

private:
  Ui::DialogVolumeFilter *ui;
  VolumeFilter* m_filter;
};

#endif // DIALOGVOLUMEFILTER_H
